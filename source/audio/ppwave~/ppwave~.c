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
	long w_begin;
	long w_len;
	float w_start;
	float w_end;
	short w_connected[2];
	t_bool w_buffer_modified;
} t_ppwave;


void *ppwave_new(t_symbol *s,  long argc, t_atom *argv);
void ppwave_free(t_ppwave *x);
t_max_err ppwave_notify(t_ppwave *x, t_symbol *s, t_symbol *msg, void *sender, void *data);
void ppwave_assist(t_ppwave *x, void *b, long m, long a, char *s);
void ppwave_limits(t_ppwave *x);
void ppwave_set(t_ppwave *x, t_symbol *s, long ac, t_atom *av);
void ppwave_float(t_ppwave *x, double f);
void ppwave_int(t_ppwave *x, long n);
void ppwave_dblclick(t_ppwave *x);
void ppwave_perform64(t_ppwave *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void ppwave_dsp64(t_ppwave *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);


static t_symbol *ps_buffer_modified;
static t_class *s_ppwave_class;


void ext_main(void *r)
{
	t_class *c = class_new("ppwave~", (method)ppwave_new, (method)ppwave_free, sizeof(t_ppwave), NULL, A_GIMME, 0);

	class_addmethod(c, (method)ppwave_dsp64,		"dsp64",	A_CANT, 0);
	class_addmethod(c, (method)ppwave_float,		"float",	A_FLOAT, 0);
	class_addmethod(c, (method)ppwave_int,		"int",		A_LONG, 0);
	class_addmethod(c, (method)ppwave_set,		"set",		A_GIMME, 0);
	class_addmethod(c, (method)ppwave_assist,		"assist",	A_CANT, 0);
	class_addmethod(c, (method)ppwave_dblclick,	"dblclick",	A_CANT, 0);
	class_addmethod(c, (method)ppwave_notify,		"notify",	A_CANT, 0);

	class_dspinit(c);
	class_register(CLASS_BOX, c);
	s_ppwave_class = c;

	ps_buffer_modified = gensym("buffer_modified");
}


void *ppwave_new(t_symbol *s,  long argc, t_atom *argv)
{
	t_ppwave *x = (t_ppwave *)object_alloc(s_ppwave_class);
	t_symbol *buf=0;
	float start=0., end=0.;
	double msr = sys_getsr() * 0.001;

	dsp_setup((t_pxobject *)x,3);
	buf = atom_getsymarg(0,argc,argv);
	start = atom_getfloatarg(1,argc,argv);
	end = atom_getfloatarg(2,argc,argv);

	x->w_name = buf;
	x->w_start = start;
	x->w_end = end;
	x->w_begin = start * msr;
	x->w_len = (end - start) * msr;
	outlet_new((t_object *)x, "signal");		// audio outlet

	// create a new buffer reference, initially referencing a buffer with the provided name
	x->w_buf = buffer_ref_new((t_object *)x, x->w_name);

	return (x);
}


void ppwave_free(t_ppwave *x)
{
	dsp_free((t_pxobject *)x);

	// must free our buffer reference when we will no longer use it
	object_free(x->w_buf);
}


// A notify method is required for our buffer reference
// This handles notifications when the buffer appears, disappears, or is modified.
t_max_err ppwave_notify(t_ppwave *x, t_symbol *s, t_symbol *msg, void *sender, void *data)
{
	if (msg == ps_buffer_modified)
		x->w_buffer_modified = true;
	return buffer_ref_notify(x->w_buf, s, msg, sender, data);
}


void ppwave_assist(t_ppwave *x, void *b, long m, long a, char *s)
{
	if (m == ASSIST_INLET) {	// inlets
		switch (a) {
		case 0:	snprintf_zero(s, 256, "(signal) Table Position (from 0 to 1)");	break;
		case 1:	snprintf_zero(s, 256, "(signal/float) Starting Table Location in ms");	break;
		case 2:	snprintf_zero(s, 256, "(signal/float) Ending Table Location in ms");	break;
		}
	}
	else {	// outlet
		snprintf_zero(s, 256, "(signal) Output %ld", a+1);
	}
}


void ppwave_limits(t_ppwave *x)
{
	t_buffer_obj *b = buffer_ref_getobject(x->w_buf); // get the actual buffer object from our reference

	if (b) {
		t_atom_long	channelcount = buffer_getchannelcount(b);		// number of floats in a frame
		t_atom_long	framecount   = buffer_getframecount(b);			// number of floats long the buffer is for a single channel
		double		msr			 = buffer_getmillisamplerate(b);	// sample rate of the buffer in samples per millisecond

		x->w_begin = (long)(x->w_start * msr) * channelcount;//buffer sr-jkc
		if (!x->w_end)	{// use entire table, eek!
			x->w_len = framecount;
		} else {
			x->w_len = (x->w_end - x->w_start) * msr; //buffer sr-jkc
		}
		// now restrict these values
		if (x->w_begin < 0)
			x->w_begin = 0;
		else if (x->w_begin >= framecount * channelcount)
			x->w_begin = (framecount - 1) * channelcount;
		if (x->w_begin + (x->w_len * channelcount) >= framecount * channelcount) {
			x->w_len = framecount - (x->w_begin / channelcount);
		}
	}
}


void ppwave_doset(t_ppwave *x, t_symbol *s, long ac, t_atom *av)
{
	t_symbol *name;
	double start, end;

	name = (ac) ? atom_getsym(av) : gensym("");
	start = (ac>1) ? atom_getfloat(av+1) : 0.;
	end = (ac>2) ? atom_getfloat(av+2) : 0.;

	if (start < 0)
		start = 0;
	if (end < 0)
		end = 0;
	x->w_start = start;
	x->w_end = end;

	buffer_ref_set(x->w_buf, name);	// change the buffer used by our buffer reference
	ppwave_limits(x);
}


// calls set the buffer ref should happen on the main thread only
void ppwave_set(t_ppwave *x, t_symbol *s, long ac, t_atom *av)
{
	defer(x, (method)ppwave_doset, s, ac, av);
}


void ppwave_float(t_ppwave *x, double f)
{
	long in = proxy_getinlet((t_object *)x);

	if (in == 1) {		// set min
		if (f < 0)
			f = 0;
		if (f > x->w_end)
			x->w_end = f;
		x->w_start = f;
		ppwave_limits(x);
	}
	else if (in == 2) {	// set max
		if (f < 0)
			f = 0;
		if (f < x->w_start)
			x->w_start = f;
		x->w_end = f;
		ppwave_limits(x);
	}
}


void ppwave_int(t_ppwave *x, long n)
{
	ppwave_float(x,(double)n);
}


void ppwave_dblclick(t_ppwave *x)
{
	buffer_view(buffer_ref_getobject(x->w_buf));
}


void ppwave_perform64(t_ppwave *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
	t_double		*in = ins[0];
	t_double		*out = outs[0];
	t_double		min = x->w_connected[0]? *ins[1] : x->w_start;
	t_double		max = x->w_connected[1]? *ins[2] : x->w_end;
	int				n = sampleframes;
	float			*b;
	long			len, dex;
	double			v;
	t_buffer_obj	*buffer = buffer_ref_getobject(x->w_buf);
	t_atom_long		channelcount;

	b = buffer_locksamples(buffer);
	if (!b)
		goto zero;

	channelcount = buffer_getchannelcount(buffer);

	if (x->w_buffer_modified) {
		x->w_buffer_modified = false;
		ppwave_limits(x);
		if (!x->w_connected[0])
			min = x->w_start;
		if (!x->w_connected[1])
			max = x->w_end;
	}

	if (min != x->w_start || max != x->w_end) {
		if (min < 0.)
			min = 0.;
		if (max < 0.)
			max = 0.;
		if (min > max)
			x->w_end = x->w_start = min;
		else {
			x->w_start = min;
			x->w_end = max;
		}
		ppwave_limits(x);
	}

	b += x->w_begin;
	len = x->w_len;

	if (channelcount == 1) {
		while (n--) {
			v = *in++;
			if (v < 0)
				v = 0;
			if (v > 1)
				v = 1;
			dex = v * (double)len;
			if (dex>len-1)
				dex = len-1;
			*out++ = b[dex];
		}
	}
	else if (channelcount > 1) {
		while (n--) {
			v = *in++;
			if (v < 0)
				v = 0;
			if (v > 1)
				v = 1;
			dex = (long)(v * (double)len) * channelcount;
			if (dex>(len-1)*channelcount)
				dex = (len-1)*channelcount;
			*out++ = b[dex];
		}
	}

	buffer_unlocksamples(buffer);
	return;
zero:
	while (n--)
		*out++ = 0.;
}


void ppwave_dsp64(t_ppwave *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
	x->w_connected[0] = count[1];
	x->w_connected[1] = count[2];
	object_method(dsp64, gensym("dsp_add64"), x, ppwave_perform64, 0, NULL);
}

