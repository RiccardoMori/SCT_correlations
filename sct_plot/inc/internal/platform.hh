#ifndef platform_h__
#define platform_h__


#define CONCATENATE_DETAIL(x, y) x##y
#define CONCATENATE(x, y) CONCATENATE_DETAIL(x, y)
#define MAKE_UNIQUE(x)  CONCATENATE(x, __LINE__)


#define  LOOP_TASK(x) auto MAKE_UNIQUE(_ret) = x; if (MAKE_UNIQUE(_ret) == return_skip) continue; else if (MAKE_UNIQUE(_ret) == return_stop) break


#define TASK_DEFINITION(TASK) template<typename next_t, typename... Args> returnTypes runTask(TASK , next_t&& next, Args&&... args)

#define LOOP_TASK_NEXT(buffer) LOOP_TASK( runTask(setBuffer(buffer, next), args...))

#define  RUN_TASK(buffer) runTask(setBuffer(buffer, next), args...)
#define RETURN_OK return return_ok
#define  RETURN_SKIP return return_skip
#define  RETURN_STOP return return_stop
#define  RETURN_NEXT_TASK(buffer) return RUN_TASK(buffer)

enum returnTypes {
  return_ok,
  return_skip,
  return_stop
};

#ifdef WIN32 
#ifndef __CINT__
#define DllExport   __declspec( dllexport )
#else
#define DllExport   
#endif // __CINT__
#else 
#define DllExport   
#endif // WIN32

#ifdef WIN32
#define SCT_FUNC __FUNCSIG__
#else
#define  SCT_FUNC __PRETTY_FUNCTION__
#endif // WIN32



#endif // platform_h__
