#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


int readfrom2DmemINT(int *x,int i,int j, int m)
{
    return x[i*m + j];
}

void deallocMem(void* array)
{
    free(array);
}


int *append(int *A, int n, int m, int *B)
{
    int j;
    int *NA;

    if(n == 0)
    {
       //NA = mallocFlat2DInt(1, 5 );

       NA = (int*) malloc((1*5)*sizeof(int)); 
    }
    else
    {
       NA = (int*) realloc(A,(n + 1)*m*sizeof(int));  
    }

    for (j=0;j < m;j++)
    {
        //printf("%d  | %d | %d \n",j,n,m);
        NA[n*m +j] = B[j];
    }
    return NA;

}

int *RegSplit(int n)
{
    int array[5] = {1,2,3,4,5};
    int *Reg,*Reg1;

    int i;
    int counts = 0;
    for(i = 0; i <n;i++)
    {
       //printf("%d \n",i);

       Reg1 = append(Reg, counts, 5, array);
       counts =  counts + 1;
       
       Reg = (int*) malloc((counts*5)*sizeof(int)); 
       memcpy(Reg, Reg1, (counts*5)*sizeof(int));
       free(Reg1);

    }

    return Reg;
    

}




int main()
{
    printf("h");
    return 1;
}
