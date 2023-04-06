/**
	@file
	ppwave~ - an oscillator using waveforms constructed by piecewise polynomials between discrete values with derivatives

	@ingroup	examples
*/

#include "ext.h"            // standard Max include, always required (except in Jitter)
#include "z_dsp.h"          // required for MSP objects
#include "math.h"
#include "ext_buffer.h"
#include "ext_atomic.h"
#include "ext_obex.h"       // required for "new" style objects


typedef struct _ppwave {
	t_pxobject w_obj;
	t_buffer_ref *w_buf;
	t_symbol *w_name;
} t_ppwave;


void *ppwave_new(t_symbol *s,  long argc, t_atom *argv);
void ppwave_free(t_ppwave *x);
void ppwave_perform64(t_ppwave *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void ppwave_dsp64(t_ppwave *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);


static t_symbol *ps_buffer_modified;
static t_class *s_ppwave_class;


void ext_main(void *r)
{
	t_class *c = class_new("ppwave~", (method)ppwave_new, (method)ppwave_free, sizeof(t_ppwave), NULL, A_GIMME, 0);

	class_addmethod(c, (method)ppwave_dsp64,		"dsp64",	A_CANT, 0);

	class_dspinit(c);
	class_register(CLASS_BOX, c);
	s_ppwave_class = c;
}


void *ppwave_new(t_symbol *s,  long argc, t_atom *argv)
{
	t_ppwave *x = (t_ppwave *)object_alloc(s_ppwave_class);

	dsp_setup((t_pxobject *)x,3);
	outlet_new((t_object *)x, "signal");		// audio outlet
    
	return (x);
}


void ppwave_free(t_ppwave *x)
{
	dsp_free((t_pxobject *)x);
}


void ppwave_perform64(t_ppwave *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
	t_double		*in = ins[0];
	t_double		*out = outs[0];
	int				n = sampleframes;

    while (n--)
    {
            *out++ = *in++;
    }
}


void ppwave_dsp64(t_ppwave *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
	object_method(dsp64, gensym("dsp_add64"), x, ppwave_perform64, 0, NULL);
}

