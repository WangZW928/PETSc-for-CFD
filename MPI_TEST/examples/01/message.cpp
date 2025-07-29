#include "mpi.h"
int main(int argc, char *argv[])
{
    char message[50];
    int messageTag = 99;
    int myrank;
    int number_amount;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if(myrank==0)
    {
        strcpy(message,"This message is from process 0");
        MPI_Send(message,strlen(message),MPI_CHAR,1,messageTag,MPI_COMM_WORLD);
        std::cout<<"message length is :"<<strlen(message)<<std::endl;
    }
    else if(myrank==1)
    {
        MPI_Status status;
        MPI_Recv(message,50,MPI_CHAR,0,messageTag,MPI_COMM_WORLD,&status);
        std::cout << "Process 1 received message: " << message << std::endl;
        MPI_Get_count(&status, MPI_CHAR, &number_amount);
        printf("rank 1 received %d numbers from 0. Message source = %d, tag = %d\n",number_amount,status.MPI_SOURCE,status.MPI_TAG);
    }
    MPI_Finalize();
    return 0;
}