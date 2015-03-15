


void complex_mult(doublecomplex* x, doublecomplex* y)
{
        doublecomplex z;
        z.re = (x->re*y->re - x->i*y->i);
        z.i = (x->i*y->re + x->re*y->i);
        x->re = z.re;
        x->i = z.i;
}
