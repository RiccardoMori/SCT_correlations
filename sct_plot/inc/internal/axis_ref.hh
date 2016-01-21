#ifndef axis_ref_h__
#define axis_ref_h__

namespace sct_corr{
  class axis_ref{
  public:
    virtual bool next() const = 0;
    virtual double get() const = 0;

  };

}
#endif // axis_h__
