#pragma once
#include <common/math.hpp>
#include <cmath>
#include <stdexcept>
namespace H02 {
template<class Scalar> Scalar viscosity(const Scalar &w,char mode) {
 using std::exp;
 if(mode!='B'&&mode!='C') throw std::invalid_argument("H02 viscosity requires B or C");
 Scalar v=w;
 if(mode=='B') { if(MPMC::scalarValue(w)<=.25) v=Scalar(0); else if(MPMC::scalarValue(w)>=.75) v=Scalar(1); else {auto s=(w-.25)/.5;v=s*s*(3.-2.*s);} }
 return exp(std::log(.002)+v*std::log(5e-5/.002));
}
}
