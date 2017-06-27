#include<Rcpp.h>
//using namespace Rcpp;

class Persona{
public:
  Persona() :edad(0), hijos(0) { };
  Persona(int edadi, int hijosi) :edad(edadi), hijos(hijosi) { };
  int edad;
  int hijos;
  void print() {
    Rcpp::Rcout << "edad = " << edad << std::endl;
    Rcpp::Rcout << "hijos = " << hijos << std::endl;
  };
  void ano(int j) {
    edad += j;
  };
};


RCPP_MODULE(Personamodule){
  Rcpp::class_<Persona>( "Persona" )
  .constructor("documentation for default constructor")
  .constructor<int,int>("documentation for constructor")
  .field( "edad", &Persona::edad, "documentation for x")
  .field( "hijos", &Persona::hijos, "documentation for y")
  .method( "print", &Persona::print, "documentation for print")
  .method( "ano", &Persona::ano, "documentation for ano");
}
