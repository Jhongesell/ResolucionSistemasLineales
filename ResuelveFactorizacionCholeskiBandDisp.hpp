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



#ifndef __ResuelveFactorizacionCholeskiBandDisp__
#define __ResuelveFactorizacionCholeskiBandDisp__


#include "ResuelveFactorizacionCholeski.hpp"
#include "MatrizBandDisp.hpp"
#include "MatrizBand.hpp"
#include "MatrizDispersa.hpp"


/// Clase para resoluci�n del sistema lineal usando factorizacion Choleski
/** @author Antonio Carrillo Ledesma
    @date primavera 2010
    @version 1.0.1
    @bug No hay errores conocidos
*/
class ResuelveFactorizacionCholeskiBandDisp: public ResuelveFactorizacionCholeski
{

private:

   /// Indica si la matriz debera ser borrada
   bool liberaMatriz;

public:


   /// Constructor de la clase
   ResuelveFactorizacionCholeskiBandDisp(void) : ResuelveFactorizacionCholeski()
   {
      M = NULL;
      RequiereMatriz = REQUIERE_MAT_DISP;
      liberaMatriz = false;
   }

   /// Constructor de la clase
   /** @param A Puntero a una matriz del tipo MatrizBandDisp */
   ResuelveFactorizacionCholeskiBandDisp(MatrizBandDisp *A) : ResuelveFactorizacionCholeski()
   {
      M = A;
      RequiereMatriz = REQUIERE_MAT_DISP;
      liberaMatriz = false;
   }

   /// Constructor de la clase
   /** @param A Puntero a una matriz del tipo MatrizBandDisp
       @param x Puntero a un Vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal */
   ResuelveFactorizacionCholeskiBandDisp(MatrizBandDisp *A, Vector *x, Vector *b) : ResuelveFactorizacionCholeski(x, b)
   {
      M = A;
      RequiereMatriz = REQUIERE_MAT_DISP;
      liberaMatriz = false;
   }

   /// Constructor de la clase
   /** @param A Puntero a una matriz del tipo MatrizBandDisp
       @param x Puntero a un Vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal */
   ResuelveFactorizacionCholeskiBandDisp(MatrizBand *A, Vector *x, Vector *b) : ResuelveFactorizacionCholeski(x, b)
   {
      M = new MatrizDispersa(A->renglones(), A->columnas(), (A->tamanoBanda() * 2 + 1));
      A->copia(M);
      RequiereMatriz = REQUIERE_MAT_BAND;
      liberaMatriz = true;
   }

   ~ResuelveFactorizacionCholeskiBandDisp()
   {
      if (liberaMatriz)
      {
         M->liberaMemoria();
         delete M;
      }
      M = NULL;
   }

   /// Factoriza la matriz A en L y U dejandolas en la misma matriz
   void factoriza(void);

   /// Resuelve el sistema lineal
   void resuelve(void);

   /// Resuelve el sistema lineal
   /** @param x Puntero a un vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal */
   void resuelve(Vector *x, Vector *b)
   {
      X = x;
      B = b;
      resuelve();
   }

};

#endif
