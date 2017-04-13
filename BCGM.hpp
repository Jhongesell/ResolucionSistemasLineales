//////////////////////////////////////////////////////////////////////////////////////////////
// Clase para resolver un sistema lineal mediante CGM                                       //
//                                                                                          //
// An�lisis y Dise�o y Programaci�n:                                                        //
//                                                                                          //
// Nombre:   Antonio Carrillo Ledesma                                                       //
// E-mail:   acl@www.mmc.geofisica.unam.mx                                                  //
// P�gina:   http://www.mmc.geofisica.unam.mx/acl                                           //
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



#ifndef __BCGM__
#define __BCGM__

#include "ResuelveSistemaLineal.hpp"
#include "MultOp.hpp"
#include "ProductoPunto.hpp"


/// Clase para resoluci�n del sistema lineal mediante CGM
/** @author Antonio Carrillo Ledesma
    @date primavera 2009
    @version 1.0.1
    @bug No hay errores conocidos
*/
class BCGM: public ResuelveSistemaLineal
{


protected:

   // Operador
   MultOp *A;

   /// Producto Punto
   ProductoPunto *prodP;

   /// N�mero m�ximo de iteraciones
   int Iter;

   /// Tolerancia
   ldouble Ep;


public:

   /// Constructor de la clase
   /** @param a Puntero a la implementaci�n a la multiplicaci�n de la matriz por el vector
       @param prod Puntero a la implementaci�n del producto punto de dos vectores
       @param iter M�ximo n�mero de interaciones
       @param ep Tolerancia m�nima */
   BCGM(MultOp &a, ProductoPunto &prod, int iter, ldouble ep) : ResuelveSistemaLineal()
   {
      A = &a;
      prodP = &prod;
      Iter = iter;
      Ep = ep;
   }

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

   /// Establece la tolerancia del metodo
   /** @param eps Tolerancia del metodo */
   void tolerancia(ldouble eps)
   {
      Ep = eps;
   }

   /// Establece el maximo numero de iteraciones
   /** @param iter Iteraciones del metodo */
   void iteraciones(int iter)
   {
      Iter = iter;
   }


};

/**
Esta clase implementa los componentes para la resoluci�n del sistema lineal mediante CGM
@example EjemploResuelveSistemaLineal.cpp */

#endif
