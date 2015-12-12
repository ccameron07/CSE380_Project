#include <masa.h>

using namespace MASA ;

int main() {
    masa_init("dumb","heateq_1d_unsteady_const") ;
    masa_disp_param() ;
    masa_init_param() ;
    
    masa_disp_param() ;

    masa_eval_1d_exact_t(1.0,1.0) ;
}
