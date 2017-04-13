//////////////////////////////////////////////////////////////////////////////////////////////
// Clase para resolver un sistema lineal usando factorizacion Choleski                      //
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



#ifndef __ResuelveFactorizacionCholeski__
#define __ResuelveFactorizacionCholeski__

#include "ResuelveSistemaLineal.hpp"
#include "MatrizBandDisp.hpp"


/// Clase para resoluci�n del sistema lineal usando factorizacion Choleski
/** @author Antonio Carrillo Ledesma
    @date primavera 2010
    @version 1.0.1
    @bug No hay errores conocidos
*/
class ResuelveFactorizacionCholeski: public ResuelveSistemaLineal
{
protected:

   /// Matriz factorizada
   bool MatrizFactorizada;

public:

   /// Constructor de la clase
   ResuelveFactorizacionCholeski(void) : ResuelveSistemaLineal()
   {
      X = NULL;
      B = NULL;
      MetodoModificaMatriz = true;
      MetodoNumerico = FACT_CHOLESKI;
      MatrizFactorizada = false;
   }

   /// Constructor de la clase
   /** @param x Puntero a un Vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal */
   ResuelveFactorizacionCholeski(Vector *x, Vector *b) : ResuelveSistemaLineal()
   {
      X = x;
      B = b;
      MetodoModificaMatriz = true;
      MetodoNumerico = FACT_CHOLESKI;
      MatrizFactorizada = false;
   }

   /// Factoriza la matriz A en L y U dejandolas en la misma matriz
   virtual void factoriza(void) = 0;

};

#endif
