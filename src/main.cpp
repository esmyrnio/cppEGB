#include "../include/keh_cst.hpp"

char eos_name[80];
double coupling, central_pressure, accuracy, relaxation;
int print_option, max_iter;


int main(int argc, char *argv[])
{
    EXECUTION::call_type = "main_call";

    for(int i=1;i<argc;i++) 
      if(argv[i][0]=='-'){
        switch(argv[i][1]){
          
          case 'f':
              sscanf(argv[i+1],"%s",eos_name);
              break;
          case 'c':
              sscanf(argv[i+1],"%lf",&coupling);
              break;
          case 'e':
              sscanf(argv[i+1],"%lf",&central_pressure);
              break;
          case 't':
              sscanf(argv[i+1],"%lf",&accuracy);
              break;
          case 'm':
              sscanf(argv[i+1],"%i",&max_iter);
              break;
          case 'l':
              sscanf(argv[i+1],"%lf",&relaxation);
              break;
          case 'p':
              sscanf(argv[i+1],"%i",&print_option);
              break;
          }
      }

    KEH_CST model(eos_name, coupling, central_pressure, accuracy, relaxation, max_iter, print_option);
    model.compute_MR();
    // std::cout<<model.mass<<" "<<model.radius<<std::endl;
    model.printModel();
    
    return 0;
}