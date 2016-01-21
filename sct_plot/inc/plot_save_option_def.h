#ifndef plot_save_option_def_h__
#define plot_save_option_def_h__


enum  plot_save_option_def
#ifndef __CINT__
  :int
#endif // !__CINT__
{
  save_to_disk,
  do_not_save_to_disk
};

#ifdef __CINT__
#pragma link C++ enum plot_save_option_def;

#endif
#endif // plot_save_option_def_h__
