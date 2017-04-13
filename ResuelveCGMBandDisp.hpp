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



#ifndef __ResuelveCGMBandada__
#define __ResuelveCGMBandada__


#include "ResuelveCGM.hpp"
#include "MatrizBandDisp.hpp"



/// Clase para resoluci�n del sistema lineal mediante CGM
/** @author Antonio Carrillo Ledesma
    @date primavera 2009
    @version 1.0.1
    @bug No hay errores conocidos
    @todo Definir clase de producto interior y poder pasarlo como argumento y user este
*/
class ResuelveCGMBandDisp: public ResuelveCGM
{

public:

   /// Constructor de la clase
   ResuelveCGMBandDisp(void) : ResuelveCGM()
   {
      M = NULL;
      C = NULL;
      RequiereMatriz = REQUIERE_MAT_BAND;
   }

   /// Constructor de la clase
   /** @param A Puntero a una matriz del tipo MatrizBandDisp */
   ResuelveCGMBandDisp(MatrizBandDisp *A) : ResuelveCGM()
   {
      M = A;
      C = NULL;
      RequiereMatriz = REQUIERE_MAT_BAND;
   }

   /// Constructor de la clase
   /** @param A Puntero a una matriz del tipo MatrizBandDisp
       @param x Puntero a un Vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal */
   ResuelveCGMBandDisp(MatrizBandDisp *A, Vector *x, Vector *b) : ResuelveCGM(x, b)
   {
      M = A;
      C = NULL;
      RequiereMatriz = REQUIERE_MAT_BAND;
   }

   /// Constructor de la clase
   /** @param A Puntero a una matriz del tipo MatrizBandDisp
       @param x Puntero a un Vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal
       @param ep Tolerancia del m�todo
       @param it N�mero m�ximo de iteraciones */
   ResuelveCGMBandDisp(MatrizBandDisp *A, Vector *x, Vector *b, ldouble ep, int it) : ResuelveCGM(x, b, ep, it)
   {
      M = A;
      C = NULL;
      RequiereMatriz = REQUIERE_MAT_BAND;
   }

};

/**
Esta clase implementa los componentes para la resoluci�n del sistema lineal mediante CGM
@example EjemploResuelveSistemaLineal.cpp */

#endif
