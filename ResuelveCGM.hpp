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


#ifndef __ResuelveCGM__
#define __ResuelveCGM__

#include "BCGM.hpp"


/// Clase para resoluci�n del sistema lineal mediante CGM standard
/** @author Antonio Carrillo Ledesma
    @date primavera 2009
    @version 1.0.1
    @bug No hay errores conocidos
*/
class ResuelveCGM: public BCGM, public  MultOp, public ProductoPunto
{

protected:

   /// Producto punto
   inline double productoPunto(Vector *u, Vector *v)
   {
      return u->productoPunto(v);
   }

   /// Multiplica Au=v
   inline void multiplica(Vector *u, Vector *v)
   {
      M->multiplica(u, v);
   }

   /// Tama�o
   inline int tamano(void)
   {
      return M->renglones();
   }


   /// Precondicionador
   Matriz *C;


public:

   /// Constructor de la clase
   ResuelveCGM(void) : BCGM(*(MultOp*) this, *(ProductoPunto*) this, 1000, 1e-5)
   {
      X = NULL;
      B = NULL;
      MetodoModificaMatriz = false;
      MetodoNumerico = CGM;
      RequiereMatriz = REQUIERE_MAT_BAND;
   }

   /// Constructor de la clase
   /** @param x Puntero a un Vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal */
   ResuelveCGM(Vector *x, Vector *b) : BCGM(*(MultOp*) this, *(ProductoPunto*) this, 1000, 1e-5)
   {
      X = x;
      B = b;
      MetodoModificaMatriz = false;
      MetodoNumerico = CGM;
      RequiereMatriz = REQUIERE_MAT_BAND;
   }

   /// Constructor de la clase
   /** @param x Puntero a un Vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal
       @param ep Tolerancia del m�todo
       @param it N�mero m�ximo de iteraciones */
   ResuelveCGM(Vector *x, Vector *b, ldouble ep, int it) : BCGM(*(MultOp*) this, *(ProductoPunto*) this, 1000, 1e-5)
   {
      X = x;
      B = b;
      MetodoModificaMatriz = false;
      MetodoNumerico = CGM;
      RequiereMatriz = REQUIERE_MAT_BAND;
   }

   /// Configura al m�todo num�rico
   /** @param ep Tolerancia del m�todo
       @param it N�mero m�ximo de iteraciones */
   inline void configuraMetodo(ldouble ep, int it)
   {
      Ep = ep;
      Iter = it;
   }
};

#endif



