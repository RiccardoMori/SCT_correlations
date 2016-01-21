#ifndef exceptions_h__
#define exceptions_h__


#include <exception>
#include <string>
#include "internal/platform.hh"

#ifndef SCT_FUNC
# define SCT_FUNC ""
#endif

#define SCT_THROWX(exc, msg) throw ::sct_corr::InitException(exc(msg), __FILE__, __LINE__, SCT_FUNC)

#define SCT_THROW(msg) SCT_THROWX(::sct_corr::Exception, (msg))


#define SCT_EXCEPTIONX(name, base) \
  class name : public base {         \
    public:                            \
   name(const std::string & msg)      \
    : base(msg) {}                   \
}

#define SCT_EXCEPTION(name) SCT_EXCEPTIONX(name, ::sct_corr::Exception)

namespace sct_corr {

class DllExport Exception : public std::exception {
public:
  Exception(const std::string & msg);
  virtual const char * what() const throw()  override {
    if (m_text.length() == 0) make_text();
    return m_text.c_str();
  }
  // This shouldn't really be const, but it must be callable on temporary objects...
  const Exception & SetLocation(const std::string & file = "",
                                unsigned line = 0,
                                const std::string & func = "") const;
  virtual ~Exception() throw() {
  }
protected:
  std::string m_msg;
private:
  void make_text() const;
  mutable std::string m_text;
  mutable std::string m_file, m_func;
  mutable unsigned m_line;
};




template <typename T>
const T & InitException(const T & e, const std::string & file, int line = 0, const std::string func = "") {
  e.SetLocation(file, line, func);
  return e;
}

// Some useful predefined exceptions
SCT_EXCEPTION(CollectionNotFound);

DllExport int handleExceptions();

}

#endif // exceptions_h__
