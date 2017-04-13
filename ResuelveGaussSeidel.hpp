//////////////////////////////////////////////////////////////////////////////////////////////
// Clase para resolver un sistema lineal mediante Gauss-Seidel                              //
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


#ifndef __ResuelveGaussSeidel__
#define __ResuelveGaussSeidel__

#include "ResuelveSistemaLineal.hpp"
#include "MatrizBandDisp.hpp"


/// Clase para resoluci�n del sistema lineal mediante Gauss-Seidel
/** @author Antonio Carrillo Ledesma
    @date primavera 2009
    @version 1.0.1
    @bug No hay errores conocidos
*/
class ResuelveGaussSeidel: public ResuelveSistemaLineal
{

protected:


   /// Tolerancia
   ldouble Ep;
   /// Iteraciones m�ximas
   int Iter;

public:


   /// Constructor de la clase
   ResuelveGaussSeidel(void) : ResuelveSistemaLineal()
   {
      X = NULL;
      B = NULL;
      Ep = 1e-5;
      Iter = 1000;
      MetodoModificaMatriz = false;
      MetodoNumerico = GAUSSSEIDEL;
   }

   /// Constructor de la clase
   /** @param x Puntero a un Vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal */
   ResuelveGaussSeidel(Vector *x, Vector *b) : ResuelveSistemaLineal()
   {
      X = x;
      B = b;
      Ep = 1e-5;
      Iter = 1000;
      MetodoModificaMatriz = false;
      MetodoNumerico = GAUSSSEIDEL;
   }

   /// Constructor de la clase
   /// Constructor de la clase
   /** @param x Puntero a un Vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal
       @param ep Tolerancia del m�todo
       @param iter N�mero m�ximo de iteraciones */
   ResuelveGaussSeidel(Vector *x, Vector *b, ldouble ep, int iter) : ResuelveSistemaLineal()
   {
      X = x;
      B = b;
      Ep = ep;
      Iter = iter;
      MetodoModificaMatriz = false;
      MetodoNumerico = GAUSSSEIDEL;
   }

};

#endif
