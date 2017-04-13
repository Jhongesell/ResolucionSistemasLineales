//////////////////////////////////////////////////////////////////////////////////////////////
// Clase para resolver un sistema lineal mediante el uso de la matriz inversa               //
// Revisi�n: 1 de Junio del 2006                                                            //
//                                                                                          //
//                                                                                          //
// An�lisis y Dise�o y Programaci�n:                                                        //
//                                                                                          //
// Nombre:   Antonio Carrillo Ledesma                                                       //
// E-mail:   acl@www.mmc.igeofcu.unam.mx                                                    //
// P�gina:   http://www.mmc.igeofcu.unam.mx/acl                                             //
//                                                                                          //
//                                                                                          //
// Este programa es software libre. Puede redistribuirlo y/o modificarlo                    //
// bajo los t�rminos de la Licencia P�blica General de GNU seg�n es                         //
// publicada por la Free Software Foundation, bien de la versi�n 2 de                       //
// dicha Licencia o bien (seg�n su elecci�n) de cualquier versi�n                           //
// posterior.                                                                               //
//                                                                                          //
// Este programa se distribuye con la esperanza de que sea �til, pero SIN                   //
// NINGUNA GARANT�A, incluso sin la garant�a MERCANTIL impl�cita o sin                      //
// garantizar la CONVENIENCIA PARA UN PROP�SITO PARTICULAR. V�ase la                        //
// Licencia P�blica General de GNU para m�s detalles.                                       //
//                                                                                          //
// Deber�a haber recibido una copia de la Licencia P�blica General junto                    //
// con este programa. Si no ha sido as�, escriba a la Free Software                         //
// Foundation, Inc., en 675 Mass Ave, Cambridge, MA 02139, EEUU.                            //
//                                                                                          //
//                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////



#ifndef __ResuelveInversa__
#define __ResuelveInversa__

#include "ResuelveSistemaLineal.hpp"
#include "MatrizDensa.hpp"


/// Clase para resoluci�n del sistema lineal mediante el uso de la matriz inversa
/** @author Antonio Carrillo Ledesma
    @date primavera 2009
    @version 1.0.1
    @bug No hay errores conocidos
*/
class ResuelveInversa: public ResuelveSistemaLineal
{
protected:

   /// Matriz factorizada
   bool MatrizInvertida;
   MatrizDensa *Inv;


public:

   /// Constructor de la clase
   ResuelveInversa(void) : ResuelveSistemaLineal()
   {
      M = NULL;
      Inv = NULL;
      X = NULL;
      B = NULL;
      MetodoModificaMatriz = true;
      MetodoNumerico = INVERSA;
      RequiereMatriz = REQUIERE_MAT_DENS;
      MatrizInvertida = false;
   }

   /// Constructor de la clase
   /** @param A Puntero a una matriz del tipo Matriz */
   ResuelveInversa(Matriz *A) : ResuelveSistemaLineal()
   {
      M = A;
      Inv = NULL;
      X = NULL;
      B = NULL;
      MetodoModificaMatriz = true;
      MetodoNumerico = INVERSA;
      RequiereMatriz = REQUIERE_MAT_DENS;
      MatrizInvertida = false;
   }

   /// Constructor de la clase
   /** @param A Puntero a una matriz del tipo Matriz
       @param x Puntero a un Vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal */
   ResuelveInversa(Matriz *A, Vector *x, Vector *b) : ResuelveSistemaLineal()
   {
      M = A;
      Inv = NULL;
      X = x;
      B = b;
      MetodoModificaMatriz = true;
      MetodoNumerico = INVERSA;
      RequiereMatriz = REQUIERE_MAT_DENS;
      MatrizInvertida = false;
   }

   /// Destructor de la clase
   ~ResuelveInversa(void)
   {
      delete Inv;
      Inv = NULL;
   }

   /// Resuelve el sistema lineal
   void resuelve(void)
   {
      if (!MatrizInvertida)
      {
         Inv = new MatrizDensa(M->renglones(), M->columnas(), "Matriz Inversa");
         invierte(M, Inv);
      }
      Inv->multiplica(B, X);
   }

   /// Resuelve el sistema lineal
   /** @param x Puntero a un vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal */
   void resuelve(Vector *x, Vector *b)
   {
      X = x;
      B = b;
      resuelve();
   }


   /// Calcula la inversa de una matriz usando el m�todo de eliminaci�n Gaussiana
   /** @param A Puntero a una matriz tipo Matriz
       @param inv Puntero a una matriz tipo Matriz que contendra la inversa */
   void invierte(Matriz *A, Matriz *inv);

};

#endif
