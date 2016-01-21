#include "internal/exceptions.hh"
#include <iostream>
#include <vector>



namespace sct_corr{

Exception::Exception(const std::string & msg)
  : m_msg(msg), m_line(0) {
}

const Exception & Exception::SetLocation(const std::string & file,
                                         unsigned line,
                                         const std::string & func) const {
  m_file = file;
  m_line = line;
  m_func = func;
  return *this;
}
std::vector<std::string> split(const std::string & str, const std::string & delim) {
  std::string s(str);
  std::vector<std::string> result;
  if (str == "") return result;
  size_t i;
  while ((i = s.find_first_of(delim)) != std::string::npos) {
    result.push_back(s.substr(0, i+1));
    s = s.substr(i + 1);
  }
  result.push_back(s);
  return result;
}
std::string text_formating(const std::string& inString) {
  auto splited = split(inString, " :/\\");
  std::string ret;
  std::string line;
  for (auto& e:splited)
  {
    
    if (line.length() + e.length()<60)
    {
      line += e;
    }else if (e.length()>=60)
    {
      line += e;
    } else {
      ret += line;
      ret += "\n  ";
      line.clear();
      line = e;
    }


   
  }
  ret += line;
  return ret;
}

void Exception::make_text() const {
  m_text = text_formating(m_msg);
  if (m_file.length() > 0) {
    m_text += "\n\n  From :\n" + text_formating(m_file);
    if (m_line > 0) {
      m_text += ":" + std::to_string(m_line);
    }
  }
  if (m_func.length() > 0) m_text += "\n\n  In: \n" + text_formating(m_func);
}

 int handleExceptions() {
  try {
    throw;
  } catch (const std::exception & e) {
    std::cerr << "Uncaught exception:\n" + std::string(e.what()) << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception:\n";
    return 1;
  }
  return 0;
}

}
