#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/utility/string_ref.hpp>
#include <boost/spirit/include/qi.hpp>
//#include <charconv>
//#include <gsl/span>
#include "EigenH5.h"
#include "highfive/highfive.hpp"
#include "blosc_filter.h"
#include "lzf/lzf_filter.h"
#include "zstd/zstd_h5plugin.h"
#include "zstd/zstd.h"
#include <H5Tpublic.h>
#include <array>
#include <stddef.h>
#include <iostream>
#include "cxxopts.hpp"








template<typename T>
class data_slice{

public:
  data_slice(const std::vector<T> & tdat, std::initializer_list<size_t>  off,std::initializer_list<size_t> chunk):
    dat(tdat),
    offset(off),
    chunksize(chunk){}
  std::vector<T> dat;
  std::vector<size_t> offset;
  std::vector<size_t> chunksize;
};


template<typename T>
class buffered_writer{
  const size_t buffer_capacity;
  std::vector<data_slice<T> > data_queue;
  HighFive::File file;
  HighFive::DataSet  ds;
public:
  buffered_writer(const std::string filename, const std::string datapath, const size_t buffer_cap):
    buffer_capacity(buffer_cap),
    file(filename,HighFive::File::ReadWrite | HighFive::File::Create),
    ds(file.getDataSet(datapath)){
    data_queue.reserve(buffer_capacity);
  }
  bool push_buffer(data_slice<T> && data_el){
    while(data_queue.size()>=buffer_capacity){
      write_data();
    }
    data_queue.emplace_back(data_el);
    return(true);
  }
  void write_data(){
    if(!data_queue.empty()){
      auto b= data_queue.back();
      ds.select(b.offset,b.chunksize,{}).write(b.dat);
      data_queue.pop_back();
    }
  }
  ~buffered_writer() {
    while(!data_queue.empty()){
      write_data();
    }
  }

};

template<typename string_type>
inline bool str_to_value(const string_type& src, double& dest)
{
    namespace qi = boost::spirit::qi;

    return qi::parse(std::cbegin(src), std::cend(src), qi::double_, dest);
}

template<typename string_type>
inline bool str_to_value(const string_type& src, int& dest)
{
    namespace qi = boost::spirit::qi;

    return qi::parse(std::cbegin(src), std::cend(src), qi::int_, dest);
}

struct data_buff{
  std::string buffer;
  std::pair<int,int> pos;
  std::string_view tbuff;
};


class mach_file{
  boost::iostreams::filtering_istream &fs;
  const std::vector<std::string> &sample_names;
  const size_t num_rows;
  const std::vector<int> &snp_indices;
  const size_t snp_ind_size;
  const size_t p;
  const size_t max_buffer_size;
  const bool SNPfirst;
  buffered_writer<double> & bw;
  //  const std::vector<size_t> line_sizes;
  boost::progress_display prog_bar;
  std::string region_buffer;
  std::string current_sample;
  //  size_t line_pos; //line cursor (should always be inside the buffer)
  size_t snp_idx; //Which snp index am I on?
  size_t line_no; // Which sample am I on?
  std::pair<int,int> buffer_pos; //Which byte in the buffer am I on (start and one past end)?
  size_t sample_offset;
  std::string	sample_id;
  std::vector<double> data;
public:
  mach_file(boost::iostreams::filtering_istream	&fs_,
	    const std::vector<std::string> &sample_names_,
	    const std::vector<int> &snp_indices_, const size_t p_,
	    const size_t buffer_size_,buffered_writer<double> & bw_,const bool SNPfirst_=true,const bool progress=true):
    fs(fs_),
    sample_names(sample_names_),
    num_rows(sample_names.size()),
    snp_indices(snp_indices_),
    snp_ind_size(snp_indices.size() ),
    p(p_),
    max_buffer_size(buffer_size_),
    bw(bw_),
    SNPfirst(SNPfirst_),
    prog_bar(num_rows)
  {
    line_no=0;
    region_buffer.reserve(max_buffer_size);
    data.reserve(max_buffer_size/6);
    buffer_pos={0,0};
    snp_idx=0;
    get_current_sample();
    //    size_t empty_rsize=name_sizes[0]+6;
  }
private:
  size_t snp_line_pos(const size_t idx) const{
    const size_t cur_snp=snp_indices[idx];
    //    const size_t empty_rsize=name_sizes[line_no]+6;
    return(cur_snp*6);
  }
  bool scan_until_snp(){
    const size_t cur_snp_pos = snp_line_pos(snp_idx);
    if(cur_snp_pos>=buffer_pos.second){
      fs.ignore(cur_snp_pos-buffer_pos.second);
      buffer_pos.second=cur_snp_pos;
      buffer_pos.first=cur_snp_pos;
      return(true);
    }else{
      size_t buff_ahead=cur_snp_pos-buffer_pos.first;
      region_buffer.substr(buff_ahead,buffer_pos.second-cur_snp_pos);
      return(false);
    }
  }

  void extend_buffer(){
    const size_t cur_snp_pos = snp_line_pos(snp_idx);
    size_t prel_snp=snp_idx;
    size_t prel_pos=snp_line_pos(prel_snp)+6;
    while((prel_snp+1)<snp_ind_size){
      if(((snp_line_pos(prel_snp+1)-cur_snp_pos)<(max_buffer_size))){
	prel_snp++;
	prel_pos=snp_line_pos(prel_snp)+6;
      }else{
	break;
      }
    }
    region_buffer.resize(prel_pos-cur_snp_pos);
      fs.read(region_buffer.data(),prel_pos-cur_snp_pos);
    buffer_pos.second=prel_pos;
  }

  void get_current_sample(){
    std::getline(fs,sample_id,'\t');
    auto ret = std::find(sample_names.begin(),sample_names.end(),sample_id);
    if(ret == sample_names.end()){
      Rcpp::stop("sample_id: "+sample_id+" not found!");
    }else{
      // std::cerr<<"sample_id: "<<sample_id<<" found on line"<<line_no<<std::endl;
      sample_offset=ret-sample_names.begin();
    }
    std::string tst;
    std::getline(fs,tst,'\t');
    // fs.read(region_buffer.data(),5);
    if(tst!="DOSE"){
      Rcpp::stop("Not at the beginning of line_no:"+std::to_string(line_no)+":\n"+tst+"\nExpecting:\nDOSE");
    }
  }

  void advance_line(){
    const size_t line_remaining=p*6-buffer_pos.second;
    fs.ignore(line_remaining);
    line_no++;
    //    region_buffer.resize(name_sizes[line_no]+6);
    get_current_sample();
    buffer_pos={0,0};
    region_buffer.clear();
    ++prog_bar;
    //  Rcpp::stop("Process Interrupted!");

    snp_idx=0;
  }
  bool parse_chunk(){
    //stuff the buffered SNPs into the data vector until:
    //1. the data vector is  full (return true)
    //2. the buffer is empty (return false)
    // If 1. We'll "artificially" advance the buffer to	the next SNP
    std::string_view tbuff(region_buffer);
    size_t buffer_size=buffer_pos.second-buffer_pos.first;
    size_t n_buffer_snps=buffer_size/6;
    size_t cur_snp_pos =snp_line_pos(snp_idx);
    size_t cur_snp=snp_indices[snp_idx];
    size_t buffer_snp=cur_snp;
    size_t start_snp=cur_snp;
    //first snp in the buffer is always cur_snp,
    //last snp in the buffer is in snp_indices

    double tres;

    for(int i=0; i<n_buffer_snps;i++){
      buffer_snp=start_snp+i;
      if(buffer_snp == cur_snp){
	str_to_value(tbuff,tres);
	data.push_back(tres);
	buffer_size-=6;
	buffer_pos.first+=6;
	tbuff=tbuff.substr(6,buffer_size);
	if((snp_idx==(snp_ind_size-1)) || (data.capacity()==data.size())){
	  write_buffer();
	  if(snp_idx==(snp_ind_size-1)){
	    if((line_no+1)<num_rows){
		advance_line();
		return(true);
	    }else{
	      return(false);
	    }
	  }
	}
	snp_idx++;
	cur_snp=snp_indices[snp_idx];
      }else{
	buffer_size-=6;
	buffer_pos.first+=6;
	tbuff=tbuff.substr(6,buffer_size);
      }
    }
    region_buffer=tbuff;
    return(true);
  }
  //  size_t sample_offset(const std::string &sample_id){
  void write_buffer(){
    const size_t data_size=data.size();
    const size_t data_start=(snp_idx>data_size) ? snp_idx-(data_size)+1 : (data_size-snp_idx-1);
    if(SNPfirst){
      bw.push_buffer(data_slice<double>(data,{data_start,sample_offset},{data_size,1}));
    }else{
      bw.push_buffer(data_slice<double>(data,{sample_offset,data_start},{1,data_size}));
    }
    data.clear();
    //    data_start+=datasize;
  }
public:
  void process_file(const bool progress=true){
    do{
      if(scan_until_snp()){
	extend_buffer();
      }
    }while(parse_chunk());

  }
};




void mach2h5(const std::string dosagefile,
             HighFive::File file,
             const std::string datapath,std::vector<int> snp_idx,
             std::vector<std::string> names,
             const int p,
             const bool SNPfirst = true,
             const size_t buffer_size=10000,
             const bool prog = false,
             const size_t buffer_vec = 1 ){
}


template<typename T>
std::vector<T> read_vec_col(const std::string filename,const std::optional<std::string> datapath=std::nullopt){
  std::vector<T> result;
  if(datapath) {
    HighFive::File file(filename, HighFive::File::ReadOnly);
    file.getDataSet(datapath.value()).read(result);
  }else{
    Rcpp::stop("Only reading from HDF5 currently supported (datapath not specified)");
  }
  return (result);
}


template<typename T>
std::vector<T> split(const std::string &s, const char delim = ',') {
    std::stringstream ss(s);
    T item;
    std::vector<T> elems;
    while (ss<<item) {
        elems.push_back(std::move(item));
        // elems.push_back(std::move(item)); // if C++11 (based on comment from @mchiasson)
    }
    return elems;
}




size_t count_p(const std::string filename, const size_t digits=5) {
    boost::iostreams::mapped_file_source mapfile;
    mapfile.open(filename);
    boost::iostreams::stream<boost::iostreams::mapped_file_source> textstream(mapfile);
    boost::iostreams::filtering_istream fs;
    fs.push(boost::iostreams::gzip_decompressor{});
    fs.push(textstream);
    std::string ts;
    std::string dose;
    std::string tsnp;
    fs >> ts >>dose >> tsnp;
    if (dose != "DOSE") {
        std::cerr<<"ts: "<<ts<<std::endl;
        Rcpp::stop("in 'count_p': expected 'DOSE' got: " + dose);
    }
    const size_t tsnp_size = tsnp.size(); // size of one SNP (should be digits)
    if (tsnp_size != digits) {
        Rcpp::stop("first dosage (" + tsnp + ") has size " + std::to_string(tsnp_size) + ", expected: " +
                   std::to_string(digits));
    }
    std::getline(fs,ts,'\n');
    const size_t p = (ts.size() / (tsnp_size + 1)) + 1;
    return p;
}


const hid_t HighFive::Filter::gzip;
const hid_t HighFive::Filter::blosc;
const hid_t HighFive::Filter::lzf;
const hid_t HighFive::Filter::zstd;
const hid_t HighFive::Filter::no_filter;




int main (int argc, char* argv[]){
  cxxopts::Options options("mach2h5", "Convert Mach VCF files to HDF5");
  options.add_options("input")
          ("m,machfile", "mach dosage file (required)",cxxopts::value<std::string>())
          ("s,SNPlist", "file/datapath to the list of SNPs to convert (required)",cxxopts::value<std::string>())
          ("n,samplenames","file/datapath to the sample IDs (required)",cxxopts::value<std::string>())
          ("p,SNPcount","total number of SNPs (in input file)",cxxopts::value<size_t>());

  options.add_options("output")
          ("f,hdf5", "HDF5 file",cxxopts::value<std::string>())
          ("d,dataset", "output dataset name", cxxopts::value<std::string>()->default_value("dosage"))
          ("c,chunksize", "chunksize dimensions",cxxopts::value<std::string>())
          ("S,SNPfirst", "Store SNP along first dimension", cxxopts::value<bool>()->implicit_value("true"))
          ("z,filter","filter to use",cxxopts::value<std::string>()->default_value("zstd"))
          ("o,filter_opts","options passed to filter",cxxopts::value<std::string>());

  options.add_options("processing")
          ("h,help", "Print help")
          ("b,buffer_size", "read buffer size",cxxopts::value<size_t>()->default_value("100000"))
          ("v,buffer_vec", "number of buffers to read at a time",cxxopts::value<size_t>()->default_value("1"));


  auto r = register_blosc(nullptr,nullptr);
  auto nr = register_zstd();
  cxxopts::ParseResult result = options.parse(argc, argv);

  if (result.count("help"))
  {
      std::cout << options.help({"", "input","output","processing"}) << std::endl;
      exit(0);
  }

  for(std::string tstring : {"machfile","SNPlist","samplenames","hdf5"}){
      if(!result.count(tstring)){
          Rcpp::stop("Option: "+tstring+" must be specified");
      }
  }




  std::string dosagefile = result["machfile"].as<std::string>();
  std::string snplistf = result["SNPlist"].as<std::string>();
  std::string samplef = result["samplenames"].as<std::string>();
  std::vector<int> snp_idx= read_vec_col<int>(snplistf,"snp_id");
  std::vector<std::string> names = read_vec_col<std::string>(samplef,"sample_id");


  std::string h5file = result["hdf5"].as<std::string>();
  std::string dataset = result["dataset"].as<std::string>();

  using namespace HighFive;
  const std::map<std::string, hid_t> filters{{"blosc", HighFive::Filter::blosc},
                                             {"no_filter", Filter::no_filter},
                                             {"none", Filter::no_filter},
                                             {"gzip", Filter::gzip},
                                             {"deflate", Filter::gzip},
                                             {"lzf", Filter::lzf},
                                             {"zstd",Filter::zstd}};
  hid_t	filt = Filter::zstd;
  if(result.count("filter")){
      std::string tfilt = result["filter"].as<std::string>();
      auto titer = filters.find(tfilt);
      if(titer == filters.end()){
          Rcpp::stop("No registered filter for filter: "+tfilt);
      }
      filt = titer->second;
  }
  size_t p=0;
  if(result.count("SNPcount")){
      p = result["SNPcount"].as<size_t>();
  }else{
      p = count_p(dosagefile);
  }

  const size_t N = names.size();
  const size_t sub_p = snp_idx.size();


  std::vector<size_t> space_dims;
  const bool SNPfirst = result["SNPfirst"].as<bool>();
  if(SNPfirst){
      space_dims = {sub_p,N};
  }else{
      space_dims = {N,sub_p};
  }

  std::vector<unsigned int> comp_opts;
  if(result.count("filter_opts")){
      comp_opts=split<unsigned int>(result["filter_opts"].as<std::string>(),',');
  }
    {
        File file(h5file, HighFive::File::ReadWrite | HighFive::File::Create |HighFive::File::Truncate);
        DataSpace space(space_dims);
        std::vector<size_t> chunk_dimensions;
        if (result.count("chunksize")) {
            chunk_dimensions = split<size_t>(result["chunksize"].as<std::string>(), ',');
        } else {
            chunk_dimensions = Filter::guess_chunk(space_dims);
        }
        HighFive::Filter filter(chunk_dimensions, filt, comp_opts);

        auto dset = file.createDataSet(dataset, space, AtomicType<double>(), filter);
        file.flush();
    }


    const size_t buffer_size = result["buffer_size"].as<size_t>();
    const size_t buffer_vec = result["buffer_vec"].as<size_t>();

    boost::iostreams::mapped_file_source mapfile;
    mapfile.open(dosagefile);
    if(!mapfile.is_open()){
        Rcpp::stop("opening	file:"+dosagefile+"failed!");
    }else{
        //std::cerr<<"File mapped, opening stream"<<std::endl;
    }
    boost::iostreams::stream<boost::iostreams::mapped_file_source> textstream(mapfile);
    boost::iostreams::filtering_istream fs;
    fs.push(boost::iostreams::gzip_decompressor{});
    fs.push(textstream);

    buffered_writer<double> bw(h5file,dataset,buffer_vec);
    mach_file mf(fs,names,snp_idx,p,buffer_size,bw,SNPfirst,true);
    mf.process_file(true);


  return 0;
}




