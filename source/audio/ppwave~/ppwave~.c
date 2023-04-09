/**
	@file
	ppwave~ - an oscillator using waveforms constructed by piecewise polynomials between discrete values with derivatives

	@ingroup	examples
*/

#include "ext.h"            // standard Max include, always required (except in Jitter)
#include "z_dsp.h"          // required for MSP objects
#include "ext_obex.h"       // required for "new" style objects


typedef struct _ppwave {
	t_pxobject w_obj;
    
    double time[2];
    double position[2];
    double velocity[2];
    double acceleration[2];
    double jerk[2];
    double snap[2];
    double crackle[2];
} t_ppwave;


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
    class_addmethod(c, (method)ppwave_time, "float",    A_FLOAT,    0);

	class_dspinit(c);
	class_register(CLASS_BOX, c);
	s_ppwave_class = c;
}


void *ppwave_new(t_symbol *s,  long argc, t_atom *argv)
{
	t_ppwave *x = (t_ppwave *)object_alloc(s_ppwave_class);

	dsp_setup((t_pxobject *)x,1);
	outlet_new((t_object *)x, "signal");		// audio outlet
    
    x->time[0]            = 1.0;
    x->position[0]        = 0.0;
    x->velocity[0]        = 0.0;
    x->acceleration[0]    = 0.0;
    x->time[1]            = 0.5;
    x->position[1]        = 1.0;
    x->velocity[1]        = 0.0;
    x->acceleration[1]    = 0.0;

	return (x);
}


void ppwave_free(t_ppwave *x)
{
	dsp_free((t_pxobject *)x);
}


void ppwave_time(t_ppwave *x, double f)
{
    if(f < 0.001) f = 0.001;
    if(f > 0.999) f = 0.999;
    x->time[1] = f;
    double tt, p0, v0, a0, pt, vt, at;

    p0 = x->position[0];
    v0 = x->velocity[0];
    a0 = x->acceleration[0];
    tt = x->time[1] - 0.0;
    pt = x->position[1];
    vt = x->velocity[1];
    at = x->acceleration[1];
    x->jerk[0] =     -3 * ( 20*(p0-pt) + 4*(3*v0+2*vt)*tt + (3*a0-  at)*tt*tt ) / (tt*tt*tt);
    x->snap[0] =     12 * ( 30*(p0-pt) + 2*(8*v0+7*vt)*tt + (3*a0-2*at)*tt*tt ) / (tt*tt*tt*tt);
    x->crackle[0] = -60 * ( 12*(p0-pt) + 6*(  v0+  vt)*tt + (  a0-  at)*tt*tt ) / (tt*tt*tt*tt*tt);

    p0 = x->position[1];
    v0 = x->velocity[1];
    a0 = x->acceleration[1];
    tt = x->time[0] - x->time[1];
    pt = x->position[0];
    vt = x->velocity[0];
    at = x->acceleration[0];
    x->jerk[1] =     -3 * ( 20*(p0-pt) + 4*(3*v0+2*vt)*tt + (3*a0-  at)*tt*tt ) / (tt*tt*tt);
    x->snap[1] =     12 * ( 30*(p0-pt) + 2*(8*v0+7*vt)*tt + (3*a0-2*at)*tt*tt ) / (tt*tt*tt*tt);
    x->crackle[1] = -60 * ( 12*(p0-pt) + 6*(  v0+  vt)*tt + (  a0-  at)*tt*tt ) / (tt*tt*tt*tt*tt);
    
//    object_post((t_object *)x, "pos was set to %.3f", x->position);
//    object_post((t_object *)x, "vel was set to %.3f", x->velocity);
//    object_post((t_object *)x, "acc was set to %.3f", x->acceleration);
//    object_post((t_object *)x, "j0 was set to %.3f", x->jerk[0]);
//    object_post((t_object *)x, "s0 was set to %.3f", x->snap[0]);
//    object_post((t_object *)x, "c0 was set to %.3f", x->crackle[0]);
//    object_post((t_object *)x, "j1 was set to %.3f", x->jerk[1]);
//    object_post((t_object *)x, "s1 was set to %.3f", x->snap[1]);
//    object_post((t_object *)x, "c1 was set to %.3f", x->crackle[1]);
}

void ppwave_perform64(t_ppwave *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
	t_double		*in = ins[0];
	t_double		*out = outs[0];
    t_double dt, dt2, dt3, dt4, dt5;
    t_double p, v, a, j, s, c;

	int	n = sampleframes;

    while (n--)
    {
        dt  = *in++;
        
        if(0 <= dt && dt < x->time[1])
        {
            dt2 = dt*dt;
            dt3 = dt2*dt;
            dt4 = dt3*dt;
            dt5 = dt4*dt;

            p = x->position[0];
            v = x->velocity[0];
            a = x->acceleration[0];
            j = x->jerk[0];
            s = x->snap[0];
            c = x->crackle[0];
        }
        else
        {
            dt -= x->time[1];
            dt2 = dt*dt;
            dt3 = dt2*dt;
            dt4 = dt3*dt;
            dt5 = dt4*dt;

            p = x->position[1];
            v = x->velocity[1];
            a = x->acceleration[1];
            j = x->jerk[1];
            s = x->snap[1];
            c = x->crackle[1];
        }
        *out++ = p + v*dt + a/2.0*dt2 + j/6.0*dt3 + s/24.0*dt4 + c/120.0*dt5;
    }
}


void ppwave_dsp64(t_ppwave *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
	object_method(dsp64, gensym("dsp_add64"), x, ppwave_perform64, 0, NULL);
}

