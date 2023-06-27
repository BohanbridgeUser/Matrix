#ifndef _BASIC_FUN_H_
#define _BASIC_FUN_H_

template<typename T>
T CSC_MAX(T a, T b)
{
    return ((a>b)? a:b);
}
template<typename T>
T CSC_MIN(T a, T b)
{
    return ((a>b)? b:a);
}
#endif