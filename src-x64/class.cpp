#include<Rcpp.h>
//using namespace Rcpp;

class Persona{
public:

  // Atributos
  int edad;
  int hijos;

  // Constructores
  Persona() :edad(0), hijos(0) { };
  Persona(int edadi, int hijosi) :edad(edadi), hijos(hijosi) { };

  //Metodos
  void print() {
    Rcpp::Rcout << "edad = " << edad << std::endl;
    Rcpp::Rcout << "hijos = " << hijos << std::endl;
  };
  void ano(int j) {
    edad += j;
  };
};


// class Sociedad{
// public:
//
//   Persona* Personas = [];
//   int length;
//
//   // Sociedad() {
//   //   Personas = NULL,
//   //   length = 0;
//   // };
//
//   void nacimiento(Persona x) {
//
//     if (Personas == NULL) {
//
//       Personas = new Persona[length + 1];
//       length++;
//       Personas[0] = x;
//
//     } else {
//
//       Persona* nuevasPersonas = new Persona[length + 1];
//
//       for(int ii = 0; ii < length; ii++) {nuevasPersonas[ii] = Personas[ii];}
//
//       nuevasPersonas[length] = x;
//       length++;
//
//       delete Personas;
//
//       Personas = nuevasPersonas;
//
//     }
//   }
// };


// void add() {
//
//   Sociedad society;
//
//   Persona juan(10, 0);
//
//   society.nacimiento(juan);
//
// }


RCPP_MODULE(Personamodule) {
  Rcpp::class_<Persona>( "Persona" )
  .constructor("documentation for default constructor")
  .constructor<int,int>("documentation for constructor")
  .field( "edad", &Persona::edad, "documentation for x")
  .field( "hijos", &Persona::hijos, "documentation for y")
  .method( "print", &Persona::print, "documentation for print")
  .method( "ano", &Persona::ano, "documentation for ano");
}


// RCPP_MODULE(Sociedadmodule) {
//   Rcpp::class_<Sociedad>( "Sociedad" )
//   .constructor("documentation for default constructor")
//   // .constructor<Persona,int>("documentation for constructor")
//   .field( "Personas", &Sociedad::Personas, "documentation for Personas")
//   .field( "length", &Sociedad::length, "documentation for length")
//   .method( "nacimiento", &Sociedad::nacimiento, "documentation for nacimiento");
// }
