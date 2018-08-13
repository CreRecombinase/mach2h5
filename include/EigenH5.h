#pragma once
#include <cmath>
#include <complex>
#include <iostream>
#include <range/v3/core.hpp>
#include <range/v3/view.hpp>
#include <range/v3/action.hpp>

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include "highfive/highfive.hpp"
#include "blosc_filter.h"
#include "lzf/lzf_filter.h"
#include "zstd/zstd_h5plugin.h"
#include "zstd/zstd.h"
#include "boost/progress.hpp"
#include <H5Tpublic.h>

template<class T> struct always_false : std::false_type {};

#if __has_include(<filesystem>)

#include <filesystem>
namespace stdx {
  using namespace ::std;
}
#elif __has_include(<experimental/filesystem>)
#   include <experimental/filesystem>
namespace stdx {
  using namespace ::std;
  using namespace ::std::experimental;
}
#else
#   error <experimental/filesystem> and <filesystem> not found
#endif

namespace Rcpp{
    void stop(std::string ms){
      std::cerr<<ms<<std::endl;
      exit(0);
    }
}