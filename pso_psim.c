/** 
* @file pso_psim.c
* @author Andre Leopoldino
* @date 2019-01-012
* @brief Modulo com implementação do algoritimo de Otimização por Enxame de Par-
* ticulas (PSO).
*
* @note Neste código fonte é utilizado o padrão de codificação do Barr Group 
* @see https://barrgroup.com/Embedded-Systems/Books/Embedded-C-Coding-Standard
* @vesion 8
*/


#include <math.h>
#include <stdlib.h>
#include "pso_psim.h"
#include <time.h>
#include <stdio.h>
#include <string.h> // para usar strcat
#include <stdint.h> // para usar uint_8


/** @todo fixed-width integers/

//#include <time.h> // para usar time

/* @see https://www.youtube.com/watch?v=xPkRL_Gt6PI&list=PLCs1uuuARiRSOIqCaXM4BL4pCIZxapUPB&index=21*/
 
/********************************************************
*   PASSOS
*
*             1) Definição do problema: minimizar o custo sendo o custo definido como o menor valor entre a potencia maxima teorica e a potencia maxima obtida
*             2) Parâmetros do PSO
*             3) Inicialização
*             4) Loop principal
*             5) Apresentação de resultados
*

Entradas:
    tensão e corrente
Saidas:
    Posição e custo de cada particula (debug)
    Posiçao e custo da melhor particula

*******************************************************/

/** @note a declaraçao de variáveis globais é para que seja possível reinicializá-las 
durante a execução da DLL pelo PSIM.
*/
static struct empty_particle particle[NPOP_MAX];
static struct best GlobalBest;
static struct gmpp_stats gmpp;
static double soma;
static float  w = W_INIT;
#if DEBUG==1
FILE *fptr;                // Arquivo de log para as saídas das partículas
#endif
static unsigned int flg_init = 0; // flag indicando se é a primeira passagem do init
static unsigned int   ite = 0;
static float vetor_pot_medias[MAX_SIZE_VETOR];/** vetor para armazenamento de
                                                potencias médias*/


static unsigned int   cnt_janela_amostra = 0;
static unsigned int   particle_id = 0; /** contador auxiliar para controlar
                                       o loop na população*/
static float          acc_ppv = 0.0f;
static float          pot_entrada_media = 0.0f;

static int            cnt_janela_estavel = 0;
unsigned int          idx; // indice genérico de 0 a 2
static char           flg_rastrear = 1;
// variaveis de entrada
volatile float        in_ipv;
volatile float        in_vpv;
static volatile float in_c1;
static volatile float in_c2;
static volatile float in_w;
static volatile unsigned int in_npop; //tamanho da populacao de particulas
static volatile char in_rst;



// variaveis para saida
static float particle_cost = 0.0f;
static float particle_position = 0.0f;
static float particle_velocity = 0.0f;
static float particle_best_position = 0.0f;
static float particle_best_cost = 0.0f;
static float global_best_cost = 0.0f;
static float global_best_position = 0.0f;
static float pot_entrada_media_out = 0.0f;

static float acc_vpv = 0.0f;  // acumulador para tensao media
static float acc_ipv = 0.0f;  // acumulador para corrente media
static float v_med_out = 0.0f;  // tensao media
static float i_med_out = 0.0f;  // tensao media

static float pot_stddev = 100.0f;

static float pot_med_mpp = POTENCIA_MAX_ARRANJO; // potencia média máxima disponível no sistema
static float d_mpp = D_INIT; // duty cicle do mpp encontrado
static float dP_mpp; // variação entre potencia no MPP e potencia média atual
static long long ite_debug = 0; // para rastrear as epocas 
static unsigned int seed;

static unsigned int max_idx_estatisticas;
static unsigned int idx_estatisticas;
/** @brief OpenSimUser
* 
* Função necessaria para integração com PSIM. Ela é executada somente na inicia-
* lização da DLL.
*
* @return Não há retorno para esta função.  
*
* @note Esta é uma nota
*/
DLLPSIM void OpenSimUser (const char *szId, 
                          const char * szNetlist,
                          void ** ptrUserData,
                          int  nInputCount,
                          int  nOutputCount, 
                          int  *pnError, 
                          char * szErrorMsg)
{  
#if DEBUG == 1
    seed = 1401591674;
    //1535037040 // levando gbest pra VARMIN

#else
    seed = (unsigned int)time (NULL);   /** para garantir que rand() sempre sejam numeros distintos*/
#endif
    seed = 1401591674;
    seed = 0;
    srand (seed);
    ite_debug = 0;
#if DEBUG==1
    // Loggando dados
    fptr = fopen (LOG_PATH, "w");
   // errno_t err;
   // err = fopen_s (&fptr,LOG_PATH, "w");

    fprintf (fptr, LOG_HEADER);
#endif
    pot_med_mpp = POTENCIA_MAX_ARRANJO; // potencia média máxima disponível no sistema
 }

/**
* @brief RunSimUser
*
* Função necessária para integração com PSIM. Ela é executada a cada ciclo de 
* processamento do PSIM.
*
* @return *out
*
* @note Nesta função que são chamados todos os cálculos e interações planta.
*/
DLLPSIM void RunSimUser (double t,
                         double delt,
                         double *in,
                         double *out,
                         void ** ptrUserData,
                         int    *pnError,
                         char   *szErrorMsg)
{
    unsigned int resto;


    //inputs
    in_npop = (unsigned int)in[2];
    in_c1 = (float)in[3];
    in_c2 = (float)in[4];
    in_w = (float)in[5];
    in_rst = (char)in[6];

    max_idx_estatisticas = QTD_EPOCH_MEDIA * in_npop - 1;
    //inicializando as variáveis pela primeira vez ou resetando

    if (flg_init == 0 || in_rst == 1)
    {
        init_particles ();
        flg_init = 1;
    }
    
    particle_position = particle[particle_id].Position;
    particle_velocity = particle[particle_id].Velocity;




    if (cnt_janela_estavel < JANELA_ESTAVEL)
    {
        cnt_janela_estavel++;
    }
    else
    {
        if (cnt_janela_amostra < JANELA_AMOSTRA)
        {
            // acumula potencia
            // freq de 100kHz 10us
            in_vpv = (float)in[0];
            in_ipv = (float)in[1];
            
            if (in_vpv < 0)
            {
                in_vpv = 0;
            }
            acc_ppv = (in_ipv*in_vpv) + acc_ppv;
            acc_vpv = acc_vpv + in_vpv;
            acc_ipv = acc_ipv + in_ipv;
            cnt_janela_amostra++;
           // out[9] = 1;// amostrando potencia
        }
        else
        {
            // calcula media da potencia e atuliza custo
            pot_entrada_media = acc_ppv / JANELA_AMOSTRA;
            v_med_out = acc_vpv / JANELA_AMOSTRA;
            i_med_out = acc_ipv / JANELA_AMOSTRA;
            particle[particle_id].Cost       = 1.0f*(POTENCIA_MAX_ARRANJO - pot_entrada_media) / POTENCIA_MAX_ARRANJO;
            particle[particle_id].AvgPower   = acc_ppv / JANELA_AMOSTRA;
            particle[particle_id].AvgVoltage = acc_vpv / JANELA_AMOSTRA;
            particle[particle_id].AvgCurrent = acc_ipv / JANELA_AMOSTRA;

            pot_entrada_media_out = pot_entrada_media;

            resto = particle_id % max_idx_estatisticas;

            //garantir janela para estatisticas (0 a max_idx_estatisticas)
            
            vetor_pot_medias[idx_estatisticas] = pot_entrada_media_out;

            acc_ppv  = 0;
            pot_entrada_media     = 0;
            cnt_janela_estavel    = 0;
            cnt_janela_amostra    = 0;
            acc_vpv               = 0;
            acc_ipv               = 0;

            particle_id++;
            // update os bests e reset de indices
            if (particle_id >= in_npop)
            {
                // update de globais
                for (idx = 0; idx < in_npop; idx++)
                {   /* Atualizando maximo pessoal*/
                    if (particle[idx].Cost < particle[idx].Best.Cost)
                    {
                        particle[idx].Best.Cost = particle[idx].Cost;
                        particle[idx].Best.Position = particle[idx].Position;

                        /* Atualizando o global best*/
                        if (particle[idx].Best.Cost < GlobalBest.Cost)
                        {
                            GlobalBest.Position = particle[idx].Best.Position;
                            GlobalBest.Cost = particle[idx].Best.Cost;
                        }
                    }
                }

                // esses valores que devem sair do código!
                global_best_cost = GlobalBest.Cost;
                global_best_position = GlobalBest.Position;

#if DEBUG==1
                // logging
                for (idx = 0; idx < in_npop; idx++)
                {                  
                    log_to_file (t, ite_debug, idx, particle[idx].Cost,
                                 particle[idx].Position, particle[idx].Velocity,
                                 particle[idx].Best.Cost,
                                 particle[idx].Best.Position,
                                 global_best_cost,global_best_position,
                                 particle[idx].AvgPower,particle[idx].AvgVoltage,
                                 particle[idx].AvgCurrent,seed);
                    // incluir o loop de atualização aqui posteriormente (após validado)
                }
#endif
                // aqui rodam os updates para o próximo ciclo!
                // @todo: depois que rodar direitinho, colocar este trecho de código no loop acima
                // onde é feito o log.
                for (idx = 0; idx < in_npop; idx++)
                {
                    /*updade velocity*/
                    /* Necessário adequar valor da velocidade para que não extrapole [0-1]
                    utilizada a função arco tangente para normalização conforme ref:
                    */
                    particle[idx].Velocity =
                        in_w * particle[idx].Velocity +
                        in_c1 * rand_local ()*(particle[idx].Best.Position - particle[idx].Position) +
                        in_c2 * rand_local ()*(GlobalBest.Position - particle[idx].Position);

                    particle[idx].Velocity = TWO_O_PI * (float)atan (particle[idx].Velocity);
                    /*update position*/
                    particle[idx].Position = particle[idx].Position + particle[idx].Velocity;

                    //garantindo [VARMIN,VARMAX]
                    particle[idx].Position = norm (particle[idx].Position);

                    /*damping inertia coeff*/
                    in_w = in_w * WDAMP;
                }

                ite++;
                if (flg_rastrear!=0)
                {
                    ite_debug++;
                }
               
                /****************************************************************
                
                                    TRECHO BUGGADO                
                
                ****************************************************************/
                
                // se buffer estiver cheio, calcular as estatisticas
               // if (ite % max_idx_estatisticas==0){
                if (idx_estatisticas == max_idx_estatisticas)
                {
                    pot_stddev = desvio_padrao(vetor_pot_medias, max_idx_estatisticas);


                    if ((pot_stddev < STD_POT_ESTAVEL) && flg_rastrear == 1)
                    {
                        // estável parar a busca!
                  
                        flg_rastrear = 0;
                        ite          = 0;
                        d_mpp        = GlobalBest.Position;
                        pot_med_mpp  = media(vetor_pot_medias, max_idx_estatisticas);

                        gmpp.ppv_avg = particle[0].AvgPower;
                        gmpp.d       = GlobalBest.Position;
                        gmpp.ipv_avg = particle[0].AvgCurrent;
                        gmpp.vpv_avg = particle[0].AvgVoltage;
                        gmpp.epoch   = ite_debug;
                        gmpp.t       = t;

                    }
                }
                dP_mpp = (float)fabs((pot_entrada_media_out - pot_med_mpp)/pot_med_mpp);
                if ((dP_mpp > DELTA_P_MIN) && (flg_rastrear == 0))
                {
                    init_particles();
                }
                //reset               
                particle_id = 0;
                
                // incrementando 
                if (idx_estatisticas < max_idx_estatisticas)
                {
                    idx_estatisticas++;
                }
                else
                {
                    idx_estatisticas = 0;
                }

            }
        }
    }

    // saídas

    if (flg_rastrear == 0)
    {
        particle_position = d_mpp;
    }

    out[0] = particle_position;
    out[1] = particle_id;
    out[2] = global_best_cost;
    out[3] = global_best_position;
    out[4] = pot_entrada_media_out; // sempre atrasada
    out[5] = particle_cost;         // sempre atrasada tbm
    out[6] = particle_best_cost;    // sempre atrasada tbm
    out[7] = particle_best_position;// sempre atrasada tbm
    out[8] = particle_velocity;
    out[9] = flg_rastrear;
    out[10] = v_med_out;//sempre atrasada
    out[11] = i_med_out;
    out[12] = seed;
    out[13] = gmpp.d;
    out[14] = gmpp.ppv_avg;
    out[15] = gmpp.vpv_avg;
    out[16] = gmpp.ipv_avg;
    out[17] = gmpp.t;
    out[18] = gmpp.epoch;
}
 
/**
* @brief CloseSimUser
* função do PSIM que é executada ao final da simulação.
*/
DLLPSIM void CloseSimUser(const char *szId, void ** ptrUserData)
{
#if DEBUG ==1
    fclose (fptr); // fechando arquivo de log
#endif
}


/**
 * @brief retorna valor aleatório entre dois números min e max.
 * 
 * @param min limite inferior 
 * @param max limite superior
 * @return float 
 * @note  srand deve ser utilizado fora da função srand(time(NULL))
 * @see http://www.cplusplus.com/forum/beginner/167036/ 
 */
float rand_in_range (float min, float max)
{
    return min + (float)(rand () / (float)(RAND_MAX + 1) * (max - min));
}

/**
* @brief rand_local
* funcao que retorna um numero float aleatorio entre VARMIN e VARMAX
* @param(in) void
* @return float aleatorio entre VARMIN e VARMAX
*/
float rand_local (void) 
{
    return (float)rand_in_range (VARMIN, VARMAX);
}

/**
 * @brief grampeamento em VARMIN e VARMAX
 * 
 * @param in valor a ser grampeado
 * @return float 
 */
float norm (float in) 
{
    float out;
    if (in > VARMAX)
    {
        out = (float)VARMAX;
    }
    if (in < VARMIN)
    {
        out = (float)VARMIN;
    }
    else
    {
        out = in;
    }    
    return (float)out;
}

#if DEBUG==1
/**
* @brief log_to_file 
* 
* função para armazenar num arquivo texto cada iteração do cógido.
*/
void log_to_file (double t,
                 int ite_id, 
                 int particle_id, 
                 float particle_cost, 
                 float particle_position,
                 float particle_velocity,
                 float particle_best_position, 
                 float particle_best_cost,
                 float global_best_cost,
                 float global_best_position,
                 float ppv_mean,
                 float vpv_mean,
                  float ipv_mean,
                  unsigned int seed)
{
     if (fptr == NULL)
    {
        exit (1);
    }
     else
     {
         fprintf (fptr, "%f; %d;%d;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%i\n", t,
                  ite_id,
                  particle_id,
                  particle_cost,
                  particle_position,
                  particle_velocity,
                  particle_best_position,
                  particle_best_cost,
                  global_best_cost,
                  global_best_position,
                  ppv_mean,
                  vpv_mean,
                  ipv_mean,
                  seed);
     }
}
#endif

void init_particles (void)
{
    float min, max;
    /** Inicializacao das globais - problema de minizacao, INF*/
    GlobalBest.Cost = 1;/* INFINITY;  Inicializando com valor infinito uma vez que
                        se trata de um problema de minização */

                        /** Inicializando população de particulas com valores aleatórios*/
    for (unsigned int idx = 0; idx < in_npop; idx++)
    {
        // inicialização aleatória por setor
        // garantindo limite mínimo em VARMIN
        if (idx == 0)
        {
            min = VARMIN;
        }
        else
        {
            min = 1.0f*idx / in_npop;
        }

        // garantindo limite superior em VARMAX
        if (idx != (in_npop - 1))
        {
            max = min + 1.0f / in_npop;
        }
        else
        {
            max = VARMAX;
        }
        //inicializando aleatoriamente por partes
        particle[idx].Position = rand_in_range (min, max);

        particle[idx].Cost = 1.0f;
        particle[idx].Velocity = TWO_O_PI*(float)atan(rand_in_range (-VARMAX, VARMAX));

        /*melhor memoria da particula neste momento....*/
        particle[idx].Best.Position = particle[idx].Position;
        particle[idx].Best.Cost = particle[idx].Cost;

        d_mpp = D_INIT;


    }

    /*@note: inicializando vetor_pot_medias com multiplos de POTENCIA_MAX_ARRANJO
    para desvio padrao elevado e todos os demais serem grandes até completar o
    preenchimento do buffer.*/
    for (unsigned int i = 0; i < in_npop; i++)
    {
        vetor_pot_medias[i] = (float)i*POTENCIA_MAX_ARRANJO;
    }

    // reiniciando variaveis
    cnt_janela_amostra = 0;
    particle_id = 0; /** contador auxiliar para controlar
                     o loop na população*/
    acc_ppv = 0.0f;
    pot_entrada_media = 0.0f;
    cnt_janela_estavel = 0;
    flg_rastrear = 1;

    // variaveis para saida
    particle_cost = 0.0f;
    particle_position = 0.0f;
    particle_velocity = 0.0f;
    particle_best_position = 0.0f;
    particle_best_cost = 0.0f;
    global_best_cost = 0.0f;
    global_best_position = 0.0f;
    pot_entrada_media_out = 0.0f;

    acc_vpv = 0.0f;  // acumulador para tensao media
    acc_ipv = 0.0f;
    v_med_out = 0.0f;  // tensao media
    i_med_out = 0.0f;  // tensao media

    pot_stddev = 100.0f;

    //pot_med_mpp = POTENCIA_MAX_ARRANJO; // potencia média máxima disponível no sistema
    d_mpp = 0.5f; // duty cicle do mpp encontrado
    w = W_INIT;

    idx_estatisticas= 0;

    gmpp.ppv_avg = 0.0f;
    gmpp.d = 0.0f;
    gmpp.ipv_avg = 0.0f;
    gmpp.vpv_avg = 0.0f;
    gmpp.epoch = 0;
    gmpp.t = 0.0f;
}


void update_swarm (double t, int ite) {

}

/**
 * @brief função para resetar particulas para nova busca de máximo global
 * em caso de alteração na característica de potência da planta
 */
void reset_particles(void)
{
}

/**
 * @brief retorna média aritmética de um vetor de tamanho n
 * 
 * @param s vetor com valores double
 * @param n quantidade de elementos do vetor
 * @return double 
 */
float media( float s[], unsigned int n )
{
    float sum = 0.0;
    unsigned int i = 0;
    for( i = 0; i < n; i++ )
    {
        sum += s[i];
    }
    return sum / n;
}

 
/**
 * @brief retorna a variância de um vetor de tamanho n
 * 
 * @param s vetor com valores double
 * @param n quantidade de elementos do vetor
 * @return double 
 */
float variancia(float s[], unsigned int n )
{
    float sum = 0.0;
    float dev = 0.0;
    float med = media(s, n );
    unsigned int i = 0;

    for( i = 0; i < n; i++ )
    {
        dev = s[i] - med;
        sum += (dev * dev);
    }

    return (float)sum / n;
}

 
/**
 * @brief desvio padrão para um vetor
 * 
 * @param s 
 * @param n 
 * @return float 
 */
float desvio_padrao(float s[], unsigned int n)
{
    float v = variancia(s, n );
    return (float)sqrt(v);
}

/**
* @brief inverter o estado de um bit
*  uso: toggle(1) return 0; toggle(0) return 1;
* @param bit
* @return char
*/
unsigned char toggle (unsigned char bit)
{
    unsigned char out;
    if (bit == 1)
    {
        out = 0;
    }
    else
    {
        out = 1;
    }
    return (unsigned char)out;
}