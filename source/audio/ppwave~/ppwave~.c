/**
	@file
	ppwave~ - an oscillator using waveforms constructed by piecewise polynomials between discrete values with derivatives

	@ingroup	examples
*/

#include "ext.h"            // standard Max include, always required (except in Jitter)
#include "z_dsp.h"          // required for MSP objects
#include "ext_obex.h"       // required for "new" style objects

#define PPWAVE_MAX_TARGETS 4

typedef struct _ppwave {
	t_pxobject w_obj;
    double target[PPWAVE_MAX_TARGETS][7];
    double temp[PPWAVE_MAX_TARGETS][7];
} t_ppwave;

typedef enum {
    PHA = 0,
    POS,
    VEL,
    ACC,
    JRK,
    SNP,
    CRK,
    TARGET_NUM_PARAMS
} target_parameter;


void *ppwave_new(t_symbol *s,  long argc, t_atom *argv);
void ppwave_free(t_ppwave *x);
void ppwave_time(t_ppwave *x, double f);
void ppwave_perform64(t_ppwave *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void ppwave_dsp64(t_ppwave *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);


static t_class *s_ppwave_class;


void ext_main(void *r)
{
	t_class *c = class_new("ppwave~", (method)ppwave_new, (method)ppwave_free, sizeof(t_ppwave), NULL, A_GIMME, 0);

	class_addmethod(c, (method)ppwave_dsp64,	"dsp64",    A_CANT,     0);
    class_addmethod(c, (method)ppwave_time,     "float",    A_FLOAT,    0);

	class_dspinit(c);
	class_register(CLASS_BOX, c);
	s_ppwave_class = c;
}


void *ppwave_new(t_symbol *s,  long argc, t_atom *argv)
{
	t_ppwave *x = (t_ppwave *)object_alloc(s_ppwave_class);

	dsp_setup((t_pxobject *)x,1);
	outlet_new((t_object *)x, "signal");		// audio outlet
    
    for (int i=0; i<PPWAVE_MAX_TARGETS; i++)
    {
        x->target[i][PHA]       = i==0 ? 1.0 : (double)i/(double)PPWAVE_MAX_TARGETS;
        x->target[i][POS]       = sin(2.0*PI*x->target[i][PHA]);
        x->target[i][VEL]       = 0.0;
        x->target[i][ACC]       = 0.0;
        x->temp[i][PHA]         = x->target[i][PHA] ;
        x->temp[i][POS]         = x->target[i][POS];
        x->temp[i][VEL]         = x->target[i][VEL];
        x->temp[i][ACC]         = x->target[i][ACC];
    }
	return (x);
}


void ppwave_free(t_ppwave *x)
{
	dsp_free((t_pxobject *)x);
}


void calculate_jerk_snap_crackle(t_ppwave *x)
{
    for (int i=0; i<PPWAVE_MAX_TARGETS; i++)
    {
        double tt, p0, v0, a0, pt, vt, at;

        for (int i=0; i<PPWAVE_MAX_TARGETS; i++)
        {
            for (int j=0; j<JRK; j++)
            {
                double diff_sign = 1.0;
                double diff = x->temp[i][j] - x->target[i][j];
                if(diff < 0.0) diff_sign = -1.0;
                if(diff * diff_sign > x->temp[i][j] * 0.001)
                {
                    x->target[i][j] += diff * 0.01;
                    //object_post((t_object *)x, "target[%i][%i] = %0.3f", i, j, x->target[i][j]);
                }
                else{
                    x->target[i][j] = x->temp[i][j];
                }
                
            }
        }
        
        // target
        int t = i+1;
        if(t >= PPWAVE_MAX_TARGETS) t = t - PPWAVE_MAX_TARGETS;
        
        // time to target
        if(i == 0)  tt = x->target[t][PHA] - 0.0;
        else        tt = x->target[t][PHA] - x->target[i][PHA];

        p0 = x->target[i][POS];
        v0 = x->target[i][VEL];
        a0 = x->target[i][ACC];
        pt = x->target[t][POS];
        vt = x->target[t][VEL];
        at = x->target[t][ACC];
        x->target[i][JRK] =  -3 * ( 20*(p0-pt) + 4*(3*v0+2*vt)*tt + (3*a0-  at)*tt*tt ) / (tt*tt*tt);
        x->target[i][SNP] =  12 * ( 30*(p0-pt) + 2*(8*v0+7*vt)*tt + (3*a0-2*at)*tt*tt ) / (tt*tt*tt*tt);
        x->target[i][CRK] = -60 * ( 12*(p0-pt) + 6*(  v0+  vt)*tt + (  a0-  at)*tt*tt ) / (tt*tt*tt*tt*tt);


    }
}

void ppwave_time(t_ppwave *x, double f)
{
    if(f < x->target[1][PHA] + 0.001) f = x->target[1][PHA] + 0.001;
    if(f > x->target[3][PHA] - 0.001) f = x->target[3][PHA] - 0.001;
//    x->target[2][PHA] = f;
    x->temp[2][PHA] = f;

    calculate_jerk_snap_crackle(x);
    
//    for (int i=0; i<PPWAVE_MAX_TARGETS; i++)
//    {
//        object_post((t_object *)x, "[i]{t,p,v,a,j,s,c} = [%i]{%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f}",
//                    i,
//                    x->target[i][PHA],
//                    x->target[i][POS],
//                    x->target[i][VEL],
//                    x->target[i][ACC],
//                    x->target[i][JRK],
//                    x->target[i][SNP],
//                    x->target[i][CRK]
//                    );
//    }

}

void ppwave_perform64(t_ppwave *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    calculate_jerk_snap_crackle(x);

	t_double		*in = ins[0];
	t_double		*out = outs[0];
    t_double dt, dt2, dt3, dt4, dt5;
    t_double p, v, a, j, s, c;

	int	n = sampleframes;

    while (n--)
    {
        double time_now  = *in++;
        
        for (int i=0; i<PPWAVE_MAX_TARGETS; i++)
        {
            // target
            int t = i+1;
            if(t >= PPWAVE_MAX_TARGETS) t = t - PPWAVE_MAX_TARGETS;
            
            // time at index and time at target
            t_double ti;
            if(i == 0)  ti = 0.0;
            else        ti = x->target[i][PHA];
            double      tt = x->target[t][PHA];
            
            if(time_now > ti && time_now < tt)
            {
                dt  = time_now - ti;
                dt2 = dt*dt;
                dt3 = dt2*dt;
                dt4 = dt3*dt;
                dt5 = dt4*dt;

                p = x->target[i][POS];
                v = x->target[i][VEL];
                a = x->target[i][ACC];
                j = x->target[i][JRK];
                s = x->target[i][SNP];
                c = x->target[i][CRK];
            }
        }
        *out++ = p + v*dt + a/2.0*dt2 + j/6.0*dt3 + s/24.0*dt4 + c/120.0*dt5;
    }
}


void ppwave_dsp64(t_ppwave *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
	object_method(dsp64, gensym("dsp_add64"), x, ppwave_perform64, 0, NULL);
}

