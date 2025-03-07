#ifndef PROLOG_H_
#define PROLOG_H_

#include <sstream>
#include <boost/format.hpp>
#include <iostream>
#include <exception>
#include <vector>

#include <Eigen/Eigen>

using namespace std;

/*
 * Logging for PROfit contained in this quick wrapper for easy command line or globel set verbosity
 */
enum log_level_t {
    LOG_CRITICAL = 0,
    LOG_ERROR = 1,
    LOG_WARNING = 2,
    LOG_INFO = 3,
    LOG_DEBUG = 4
};

extern log_level_t GLOBAL_LEVEL;

namespace log_impl {


    class formatted_log_t {
        public:
            formatted_log_t( log_level_t level, const wchar_t* msg ) : level(level), fmt(msg) {}
            ~formatted_log_t() {
                // GLOBAL_LEVEL is a global variable and could be changed at runtime
                // Any customization could be here
                if ( level <= GLOBAL_LEVEL ) wcout << level << L" " << fmt << endl;
            }        
            template <typename T> 
                formatted_log_t& operator %(T value) {
                    fmt % value;
                    return *this;
                }
            template <typename T>
            formatted_log_t& operator %(const std::vector<T>& vec) {
                std::wstringstream ss;
                ss << L"[";
                for (size_t i = 0; i < vec.size(); ++i) {
                    if (i != 0) ss << L", ";
                    ss << vec[i];
                }
                ss << L"]";
                fmt % ss.str();
                return *this;
            }
            template<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime, int Options,
                        int MaxRowsAtCompileTime, int MaxColsAtCompileTime>
            formatted_log_t& operator %(const Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>& vec) {
                std::wstringstream ss;
                ss << L"[";
                // TODO: Probably broken for matrices
                for (int i = 0; i < vec.size(); ++i) {
                    if (i != 0) ss << L", ";
                    ss << vec[i];
                }
                ss << L"]";
                fmt % ss.str();
                return *this;
            }

        protected:
            log_level_t     level;
            boost::wformat      fmt;
    };

    template <>
    inline formatted_log_t& formatted_log_t::operator %(const std::vector<std::string>& vec) {
        std::wstringstream ss;
        ss << L"[";
        for (size_t i = 0; i < vec.size(); ++i) {
            if (i != 0) ss << L", ";
            ss << vec[i].c_str();
        }
        ss << L"]";
        fmt % ss.str();
        return *this;
    }
}//namespace log_impl
// Helper function. Class formatted_log_t will not be used directly.
template <log_level_t level>
log_impl::formatted_log_t log(const wchar_t* msg) {
    return log_impl::formatted_log_t( level, msg );
}

#endif
