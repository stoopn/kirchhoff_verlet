#include "wire.h"
//#include <alps/osiris/comm.h>

int main(int argc, char** argv) {

if (argc<2) {
	std::cout<<"Please supply an input file"<<std::endl;
	exit(0);
}
std::cout<<"Starting mywsim with "<<std::string(argv[1])<<"\n";
    std::string filen=std::string(argv[1]);
    WireSim mywsim(filen);
    
    bool do_run=true;
    int checkblock=100000;
    int ctr=0;
    while (do_run)
    {
        mywsim.dostep();
        ctr++;
        if (ctr>checkblock) {
            if (mywsim.work_done()==1.0) do_run=false;
            ctr=0;
            std::cout<<"Checking if finished..\n";
        }
        
    }
std::cout<<"En d\n";

}
