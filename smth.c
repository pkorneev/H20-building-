#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <semaphore.h>

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/wait.h>


typedef struct {
    int no;                // number of oxygens
    int nh;                // number of hydrogens
    int ti;                // timeout for atom whaiting
    int tb;                // timeout for bonding
    int action_counter;    // action counter for print in output
    int mol_id;            // id of molecule
    int oxygens;           // amount of oxygens
    int hydrogens;         // amount of hydrogens
    int mols;
    int h_cnt;
    int o_cnt;
   
    

    sem_t o_queue_sem;     // for queue oxyg
    sem_t h_queue_sem;     // for queue hydro
    sem_t print_sem;       // for printout
    sem_t oxygens_waiting; // wait
    sem_t hydrogens_waiting;
    sem_t thats_all_sem;   // end 
    sem_t no_hydrogens;    // no remaining  
    sem_t no_oxygens;      //
    sem_t oxygen_create;   // create molecule oxygen
    sem_t hydrogen_create; // create molecule hydrogen
    sem_t h_cnt_sem;       // cnt created oxygens;
    sem_t o_cnt_sem;       // cnt created hydrogens;
    


} Shared_mem;

int shm_id = 0; // global values only for shmem
Shared_mem *shm; // our shmem

void error(char *error, int err_code){
    fprintf(stderr, "%s\n", error);
    exit(err_code);
}
// parse cl arguments     NO: oxygen quantity
//                        NH: hydrogen quantity
void parse_args(int argc, char*argv[]) 
{
   if(argc != 5){
       fprintf(stderr,"ERROR!\nSyntax of start: ./proj2 NO NH TI TB\n"); // if we didn't get 4 arguments to our program -> exit  
       exit(EXIT_FAILURE);
   }
   char *endptr;
   shm->no = strtol(argv[1], &endptr, 10);
   if(endptr == argv[1] || *endptr != '\0')     error("ERROR!\nwrong args", EXIT_FAILURE);
   shm->nh = strtol(argv[2], &endptr, 10);
   if(endptr == argv[2] || *endptr != '\0')     error("ERROR!\nwrong args", EXIT_FAILURE);
   shm->ti = strtol(argv[3], &endptr, 10);
   if(endptr == argv[3] || *endptr != '\0')     error("ERROR!\nwrong args", EXIT_FAILURE);
   shm->tb = strtol(argv[4], &endptr, 10);
   if(endptr == argv[4] || *endptr != '\0')     error("ERROR!\nwrong args", EXIT_FAILURE);

   shm->no = atoi(argv[1]);
   shm->nh = atoi(argv[2]);
   shm->ti = atoi(argv[3]);
   shm->tb = atoi(argv[4]);
   

   if(shm->ti < 0 || shm->ti > 1000 || shm->tb < 0 || shm->tb > 1000){
       fprintf(stderr,"ERROR!\nSyntax of start: ./proj2 NO NH TI TB\nTI is 0<=TI<=1000\nTB is 0<=TB<=1000\n");
       exit(EXIT_FAILURE);
   }
   if(shm->no < 0 || shm->nh < 0){
       fprintf(stderr,"ERROR!\nNO and NH must be > 1");
       exit(EXIT_FAILURE);
   }
   return;
}
//create shmem for our processes (Shared_mem* shm)
void create_shared_memory()
{
    
    shm_id = shmget(IPC_PRIVATE, sizeof(Shared_mem), IPC_CREAT | 0666);
    if(shm_id < 0){
        fprintf(stderr,"ERROR!\nShared memory couldn't be created.\n");
        exit(EXIT_FAILURE);
    }
    shm = (Shared_mem*)shmat(shm_id, NULL, 0);
    if(shm == (Shared_mem *) -1){
        fprintf(stderr,"ERROR!\nShared memory couldn't be created.\n");
        exit(EXIT_FAILURE);
    }
    return;
}



void delete_shared_memory () {  // destroy semaphores
    if (sem_destroy(&shm->print_sem) == -1)              error("ERROR!\nsem_destroy failed",EXIT_FAILURE);
    if (sem_destroy(&shm->o_queue_sem) == -1)            error("ERROR!\nsem_destroy failed",EXIT_FAILURE);
    if (sem_destroy(&shm->h_queue_sem) == -1)            error("ERROR!\nsem_destroy failed",EXIT_FAILURE);
    if (sem_destroy(&shm->thats_all_sem) == -1)          error("ERROR!\nsem_destroy failed",EXIT_FAILURE);
    if (sem_destroy(&shm->no_hydrogens) == -1)           error("ERROR!\nsem_destroy failed",EXIT_FAILURE);
    if (sem_destroy(&shm->no_oxygens) == -1)             error("ERROR!\nsem_destroy failed",EXIT_FAILURE);
    if (sem_destroy(&shm->oxygen_create) == -1)          error("ERROR!\nsem_destroy failed",EXIT_FAILURE);
    if (sem_destroy(&shm->hydrogen_create) == -1)        error("ERROR!\nsem_destroy failed",EXIT_FAILURE);
    if (sem_destroy(&shm->hydrogens_waiting) == -1)      error("ERROR!\nsem_destroy failed",EXIT_FAILURE);
    if (sem_destroy(&shm->oxygens_waiting) == -1)        error("ERROR!\nsem_destroy failed",EXIT_FAILURE);
    if (sem_destroy(&shm->o_cnt_sem) == -1)              error("ERROR!\nsem_destroy failed",EXIT_FAILURE);
    if (sem_destroy(&shm->h_cnt_sem) == -1)              error("ERROR!\nsem_destroy failed",EXIT_FAILURE);

    if(shmdt(shm)==-1)
    {
        fprintf(stderr,"err_print!\nshmdt failed.\n");
        return;
    }
    shmctl(shm_id, IPC_RMID, NULL);
}

FILE* open_file() // open file
{
    FILE *f = fopen("proj2.out", "w");
    if(f == NULL)
    {
        fprintf(stderr, "ERROR!\nProblem with openning file\n");
        exit(EXIT_FAILURE);
    }
    return f;
}

void sleeping(int arg) // sleeping for arg msec. ( TI or TB)
{
    int timeout = (rand() % (arg + 1)) * 1000; 
    usleep(timeout);
}

void sem_initialize(){  // initialize semaphores 
    if (sem_init(&shm->print_sem,1,1) == -1)                  error("ERROR!\nsem_init failed",EXIT_FAILURE);
    if (sem_init(&shm->h_queue_sem,1,0) == -1)                error("ERROR!\nsem_init failed",EXIT_FAILURE);
    if (sem_init(&shm->o_queue_sem,1,1) == -1)                error("ERROR!\nsem_init failed",EXIT_FAILURE);
    if (sem_init(&shm->hydrogens_waiting,1,0) == -1)          error("ERROR!\nsem_init failed",EXIT_FAILURE);
    if (sem_init(&shm->oxygens_waiting,1,0) == -1)            error("ERROR!\nsem_init failed",EXIT_FAILURE); 
    if (sem_init(&shm->thats_all_sem,1,0) == -1)              error("ERROR!\nsem_init failed",EXIT_FAILURE);
    if (sem_init(&shm->no_hydrogens,1,0) == -1)               error("ERROR!\nsem_init failed",EXIT_FAILURE);
    if (sem_init(&shm->no_oxygens,1,0) == -1)                 error("ERROR!\nsem_init failed",EXIT_FAILURE);
    if (sem_init(&shm->hydrogen_create,1,0) == -1)            error("ERROR!\nsem_init failed",EXIT_FAILURE);
    if (sem_init(&shm->oxygen_create,1,0) == -1)              error("ERROR!\nsem_init failed",EXIT_FAILURE);
    if (sem_init(&shm->o_cnt_sem,1,1) == -1)                  error("ERROR!\nsem_init failed",EXIT_FAILURE);
    if (sem_init(&shm->h_cnt_sem,1,1) == -1)                  error("ERROR!\nsem_init failed",EXIT_FAILURE);
    return;
}

void set_value(){ //setup nedeed values
    shm->action_counter = 1;
    shm->mol_id = 0;
    shm->hydrogens = shm->nh;
    shm->oxygens = shm->no;
    shm->mols = 0;
    shm->o_cnt = 0;
    shm->h_cnt = 0;

}

void print(FILE *f, char* s, int n1, int n2) { //print with mutexes and inc actioncounter
    sem_wait(&shm->print_sem);
    fprintf(f,"%d: ", shm->action_counter);
    fprintf(f, s, n1, n2);
    shm->action_counter++;
    sem_post(&shm->print_sem);	
}

void posting(){ // post everything, where somebody could somehow be stopped 
    sem_post(&shm->o_queue_sem);
    sem_post(&shm->h_queue_sem);
    sem_post(&shm->oxygens_waiting);
    sem_post(&shm->no_hydrogens);
    sem_post(&shm->no_oxygens);
    sem_post(&shm->thats_all_sem);
}

void creating(int flag, int id, FILE *f){ // flag 1-for oxygen, 2-for hydrogen
    if (flag == 1){
    shm->mol_id += 1;
    sem_wait(&shm->hydrogens_waiting);
    sem_wait(&shm->hydrogens_waiting);
    sem_post(&shm->oxygens_waiting);
    sem_post(&shm->oxygens_waiting);
    print(f,"O %d: creating molecule %d\n", id, shm->mol_id);
    sem_post(&shm->h_queue_sem);
    sem_post(&shm->h_queue_sem);
    sem_wait(&shm->hydrogen_create);
    sem_wait(&shm->hydrogen_create);
    sleeping(shm->tb);
    print(f,"O %d: molecule %d created\n", id, shm->mol_id);

    sem_post(&shm->oxygen_create);
    sem_post(&shm->oxygen_create);
    sem_wait(&shm->hydrogen_create);
    sem_wait(&shm->hydrogen_create);
    
    shm->h_cnt -= 2;
    shm->o_cnt -= 1;
    shm->hydrogens -= 2;  
    shm->oxygens -= 1;

    sem_post(&shm->o_queue_sem);
    } else if (flag == 2){
    print(f,"H %d: creating molecule %d\n", id, shm->mol_id);
    sem_post(&shm->hydrogen_create);
    sem_wait(&shm->oxygen_create);
    print(f,"H %d: molecule %d created\n", id, shm->mol_id);
    }
}

void process_oxygen(FILE *f, int o_id){
    srand(getpid() * time(NULL));
    sem_wait(&shm->o_cnt_sem);
    shm->o_cnt += 1;
    sem_post(&shm->o_cnt_sem);
    print(f, "O %d: started\n", o_id, 0);
    sleeping(shm->ti);
    print(f, "O %d: going to queue\n", o_id, 0);
    if(o_id == shm->no) sem_post(&shm->no_oxygens);
    sem_wait(&shm->o_queue_sem);
    if(shm->hydrogens < 2){ // when not enough hydrogens to bonding
        sem_wait(&shm->thats_all_sem);
        sem_wait(&shm->no_hydrogens);
        sem_wait(&shm->no_oxygens);
        print(f,"O %d: not enough H\n", o_id, 0);
        posting();
        exit(0);
    }
    creating(1, o_id, f);
    if (shm->oxygens < 1)
    {
        sem_post(&shm->h_queue_sem);
        sem_post(&shm->oxygens_waiting);
        sem_post(&shm->thats_all_sem);
    } else if(shm->oxygens >= 1 && shm->hydrogens < 2) sem_post(&shm->thats_all_sem);

    exit(0);
            


}

void process_hydrogen(FILE *f, int h_id){
    srand(getpid() * time(NULL));
    sem_wait(&shm->h_cnt_sem);
    shm->h_cnt += 2;
    sem_post(&shm->h_cnt_sem);
    print(f, "H %d: started\n", h_id, 0);
    sleeping(shm->ti);
    print(f, "H %d: going to queue\n", h_id, 0);
    if(h_id == shm->nh) sem_post(&shm->no_hydrogens);
    sem_post(&shm->hydrogens_waiting);
    sem_wait(&shm->oxygens_waiting);
    sem_wait(&shm->h_queue_sem);
    if(shm->hydrogens < 2 || shm->oxygens < 1){ // when not enough atoms to bonding
        sem_wait(&shm->thats_all_sem);
        sem_wait(&shm->no_hydrogens);
        sem_wait(&shm->no_oxygens);
        print(f,"H %d: not enough O or H\n", h_id, 0);
        posting();
        exit(0);
    }
    creating(2, h_id, f);
    if (shm->hydrogens < 2)
    {
        sem_post(&shm->thats_all_sem);
    }
    sem_post(&shm->hydrogen_create);
    
    exit(0);
           
}


int main(int argc, char *argv[])
{
    create_shared_memory();
    parse_args(argc, argv);
    FILE *f = open_file();
    setbuf(f, 0); 

    set_value();
    sem_initialize();  
    
    pid_t pidsO[shm->no];
    pid_t pidsH[shm->nh];
    if(shm->no == 0 && shm->nh == 0)error("ERROR!\n0 atoms\nwrong args\n",EXIT_FAILURE);
    if(shm->no == 0){
        sem_post(&shm->no_oxygens);
        sem_post(&shm->h_queue_sem);
        //sem_post(&shm->o_queue_sem);
        sem_post(&shm->oxygens_waiting);
        sem_post(&shm->thats_all_sem);
    }
    if(shm->nh < 2){
        sem_post(&shm->no_hydrogens);        
        sem_post(&shm->h_queue_sem);
        sem_post(&shm->oxygens_waiting);
        sem_post(&shm->thats_all_sem);
    } 

    for (int i = 0; i < shm->no; i++){
        pidsO[i] = fork();
        if(pidsO[i] == 0){
            process_oxygen(f, (i + 1));
        } else if(pidsO[i] == -1){
            error("ERROR!\nfork failed\n",EXIT_FAILURE);
        }  
    }

    for (int i = 0; i < shm->nh; i++){
        pidsH[i] = fork();
        if(pidsH[i] == 0){
            process_hydrogen(f, (i + 1));
        } else if(pidsH[i] == -1){
            error("ERROR!\nfork failed\n",EXIT_FAILURE);
        }
    }
    

    while(wait(NULL) > 0);        
    delete_shared_memory();
    fclose(f);
    exit(0);

    return 0;
}