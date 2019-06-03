/** 
* @file pso_psim.h
* @author Andre Leopoldino
* @date 2018-03-06
* @brief Definições de constantes e esqueletos de funções obrgatorias para o PSIM
*/
#pragma warning(suppress : 4996)

#define DEBUG 1

#ifndef BUILD_DLL_PSIM
    #define DLLPSIM __declspec(dllexport)
#else
    #define DLLPSIM __declspec(dllimport)
#endif

/* Funções principais*/

/**
 * @brief OpenSimUser
 * 
 * @param szId 
 * @param szNetlist 
 * @param nInputCount 
 * @param nOutputCount 
 * @param pnError 
 * @param szErrorMsg 
 * @return DLLPSIM OpenSimUser 
 */
DLLPSIM void OpenSimUser (const char *szId,
                          const char * szNetlist,
                          void ** ptrUserData,
                          int  nInputCount,
                          int  nOutputCount,
                          int  *pnError,
                          char * szErrorMsg);

/**
 * @brief 
 * 
 * @param t 
 * @param delt 
 * @param in 
 * @param out 
 * @param pnError 
 * @param szErrorMsg 
 * @return DLLPSIM RunSimUser 
 */
DLLPSIM void RunSimUser (double t,
                         double delt,
                         double *in,
                         double *out,
                         void ** ptrUserData,
                         int    *pnError,
                         char   *szErrorMsg);
DLLPSIM void CloseSimUser (const char *szId, void ** ptrUserData);

/**
* @brief best
*
* Prototipo da estrutura para armazenar os melhores valores.
*/
struct best
{
    float Cost;      /** Melhor custo em tempo de execução da partícula */
    float Position;  /** Melhor posição em tempo de execução da partícula */
} best;

/** 
* @brief empty_particle
* 
* prototipo da estrutura da partícula com as características necessárias
*/
struct empty_particle
{
    float Position; /** Posição atual da partícula */
    float Velocity; /** Velocidade atual da partícula */
    float Cost;     /** Velocidade atual da partícula */
    float AvgPower;
    float AvgVoltage;
    float AvgCurrent;
    /* estrutura dentro da estrutura*/
    struct best Best;
};

struct gmpp_stats
{
    float d;
    float ppv_avg;
    float vpv_avg;
    float ipv_avg;
    float t;
    float epoch;
};


/* Parametros */
#define MAXIT          (int)20	      /** maximo de epocas */
#define WDAMP          (float)1       // dumping ratio of inertia coeff
#define NPOP_MAX       (unsigned int)11
#define VARMIN         (float)0.001  // duty cicle minima de busca
#define VARMAX         (float)0.999  // duty cicle maximo de busca
#define JANELA_AMOSTRA (int)50       /** tamanho da janela para acumular os 
                                        resultados e tirar medias 50*/
#define LOG_PATH       "C:\\pso_psim_dll.csv"
#define LOG_HEADER     "exec_time; ite_id; particle_id; particle_cost; particle_position; particle_velocity; particle_best_cost;particle_best_position; global_best_cost; global_best_position;ppv_media;vpv_media;ipv_media;seed\n"
#define POTENCIA_MAX_ARRANJO (float)400.0// potência máxima do arranjo de paineis
#define THERESHOLD_UPDATE (float)0.01  /* para que seja realizada uma atualização 
                                        o valor novo deve ser 5%inferior ao anterior*/
#define JANELA_ESTAVEL (int)450       /** quantidade de ciclos para estabilizar 
                                        a atuação do dutycicle*/
#define PI             (float)3.14159265358979323846
#define TWO_O_PI       (float)2.0/PI
#define DELTA_P_MIN    (float)0.05
#define MAX_SIZE_VETOR (unsigned int)3*NPOP_MAX       // Tamanho max do vetor de acumulo de 
#define QTD_EPOCH_MEDIA (unsigned int) 1

#define STD_POT_ESTAVEL (float)0.8// identificando estabilidade na potencia
#define W_INIT         (float)0.729
#define D_INIT         (float)0.7
/*******************************************************************************
 *                        prototipos das funções
 ******************************************************************************/

float rand_in_range (float, float);
float rand_local (void);
float norm (float);
void  init_particles (void); // resetar variaveis aleatorias
void  update_swarm (double, int);
void  reset_particles(void);
float media(float s[], unsigned int  );
float variancia(float s[], unsigned int);
float desvio_padrao(float s[], unsigned int);
unsigned char toggle (unsigned char);

#if DEBUG==1
    void  log_to_file (double, int, int, float, float, float, float, float, float, float, float, float, float, unsigned int);
#endif
