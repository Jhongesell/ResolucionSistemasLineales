//////////////////////////////////////////////////////////////////////////////////////////////
// Clase para resolver un sistema Ax=b                                                      //
// Revisi�n: 20 de Junio del 2006                                                           //
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

#include <math.h>
#include "ResuelveGaussSeidelBandDisp.hpp"


// Resuelve Ax=b usando el metodo Gauss-Seidel
void ResuelveGaussSeidelBandDisp::resuelve(void)
{
#ifdef DEPURAR
   if (!M->matrizCuadrada()) error("matriz no cuadrada");
#endif

   int i, j, ind, xcol;
   ldouble sum, mx1, m1;
   Vector *xt = new Vector(M->renglones());
   if (!xt) error("no hay memoria para el vector auxiliar");

   for (NumIteraciones = 1; NumIteraciones <= Iter; NumIteraciones++)
   {
      for (i = 0; i < M->renglones(); i++)
      {
         if (M->retorna(i, i) == 0.0)
         {
            delete xt;
            return;
         }
         sum = 0.0;
         xcol = M->retornaNumeroColumnasBanda(i);
         for (ind = 0; ind < xcol; ind ++)
         {
            j = M->retornaNumeroColumna(i, ind);
            if (i == j) continue;
            sum += M->retorna(i, j) * xt->retorna(j);
            //~ if (j < i) sum += M->retorna(i,j) * xt->retorna(j);
            //~ if (j >= i+1) sum += M->retorna(i,j) * X->retorna(j);
         }
         xt->asigna(i, (B->retorna(i) - sum) / M->retorna(i, i));
      }
      mx1 = 0.0;
      for (i = 0; i < M->renglones(); i++)
      {
#ifdef __Double__
         m1 = fabs(xt->retorna(i) - X->retorna(i));
#else
         m1 = fabsl(xt->retorna(i) - X->retorna(i));
#endif
         if (m1 > mx1) mx1 = m1;
         X->asigna(i, xt->retorna(i));
      }
      //~ printf("\n%f",mx1);
      if (mx1 <= Ep) break;
   }
   delete xt;
}
