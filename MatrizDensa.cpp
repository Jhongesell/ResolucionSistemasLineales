//////////////////////////////////////////////////////////////////////////////////////////////
// Clase para el trabajar con matrices densas de punto flotante                             //
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

#include "MatrizDensa.hpp"


// Solicita la memoria necesaria para contener los valores de la matriz
// @param ren N�mero de renglones de la matriz
// @param col N�mero de columnas de la matriz
void MatrizDensa::solicitaMemoria(const int ren, const int col)
{
   int i, j;
   M = new ldouble *[ren];
   if (M == NULL) faltaMemoria();

   for (i = 0; i < ren; i++)
   {
      M[i] = new ldouble[col];
      if (M[i] == NULL) faltaMemoria();
   }

   Col = col;
   Ren = ren;
   Tipo_Matriz = MATRIZ_DENSA;

   for (i = 0; i < Ren; i++)
   {
      for (j = 0; j < Col; j++) M[i][j] = 0;
   }

}


// Libera la memoria solicitada para la matriz
void MatrizDensa::liberaMemoria(void)
{
   int i;
   if (M)
   {
      for (i = 0; i < Ren; i++)
      {
         delete []M[i];
         M[i] = NULL;
      }
      delete []M;
      M = NULL;
   }
}


// Inicializa la matriz al valor indicado
// @param val Valor por omisi�n para inicializar la matriz
void MatrizDensa::inicializa(const ldouble val)
{
   int i, j;
   for (i = 0; i < Ren; i++)
   {
      for (j = 0; j < Col; j++) M[i][j] = val;
   }
}



#ifdef DEPURAR

// Asigna el valor indicado en el renglo y columna solicitado
// @param ren Renglon
// @param col Columna
// @param val Valor
void MatrizDensa::asigna(const int ren, const int col, const ldouble val)
{
   if (ren < 0 || ren >= Ren || col < 0 || col >= Col)
   {
      printf("Error al asignar, indices desbordados (%d, %d)\n", ren, col);
      visualiza(0);
      exit(1);
   }
   M[ren][col] = val;
}


// Retorna el valor del renglon y columna solicitado
// @param ren Renglon
// @param col Columna
// @return Valor
ldouble MatrizDensa::retorna(const int ren, const int col)
{
   if (ren < 0 || ren >= Ren || col < 0 || col >= Col)
   {
      printf("Error al recuperar, indices desbordados (%d, %d)\n", ren, col);
      visualiza(0);
      exit(1);
   }

   return M[ren][col];
}

#endif



// Multiplica la matriz por el escalar pasado como parametro
// @param esc Escalar
void MatrizDensa::multiplica(ldouble esc)
{
   int i, j;

   for (i = 0; i < Ren; i++)
   {
      for (j = 0; j < Col; j++) asigna(i, j, retorna(i, j) * esc);
   }
}


// Multiplica la matriz A por la matriz B
// @param a Puntero a matriz densa
// @param b Puntero a matriz densa
void MatrizDensa::multiplica(MatrizDensa *a, MatrizDensa *b)
{
#ifdef DEPURAR
   if (a->columnas() != b->renglones())
   {
      printf("Error en dimensi�n de las matrices para su multiplicarci�n\n");
      visualizaInformacion();
      a->visualizaInformacion();
      b->visualizaInformacion();
      exit(1);
   }
#endif

   int i, j, k;
   ldouble v;
   for (i = 0; i < a->renglones(); i++)
   {
      for (j = 0; j < b->columnas(); j++)
      {
         v = 0;
         for (k = 0; k < a->columnas(); k++)
         {
            v = v + a->retorna(i, k) * b->retorna(k, j);
         }
         asigna(i, j, v);
      }
   }
}


// Multiplica la matriz por el vector B dejando el Resultado en R
// @param b Puntero a un Vector
// @param r Puntero a un Vector
void MatrizDensa::multiplica(Vector *b, Vector *r)
{
#ifdef DEPURAR
   if (columnas() != b->columnas() || renglones() != r->columnas())
   {
      printf("Multiplica Mb=r Error en dimensi�n\n");
      visualizaInformacion();
      b->visualizaInformacion();
      r->visualizaInformacion();
      exit(1);
   }
#endif


   int i, k;
   ldouble v;


   for (i = 0; i < renglones(); i++)
   {
      v = 0;
      for (k = 0; k < columnas(); k++)
      {
         v = v + retorna(i, k) * b->retorna(k);
      }
      r->asigna(i, v);
   }
}




