#include "Random.h"
#include <math.h>

using namespace std;

double Random::randn ()
{
//   static bool _havesmpl = false;
//   static double _smpl;

  if(_havesmpl)
  {
	  _havesmpl = false;
	  return _smpl;
  }
  else
  {
	  double x1, x2, w;

	  do
		{
		  x1 = 2.0 * (rand_r (&_seed) / (double) RAND_MAX) - 1.0;
		  x2 = 2.0 * (rand_r (&_seed) / (double) RAND_MAX) - 1.0;
		  w = x1 * x1 + x2 * x2;
		}
	  while (w >= 1.0);

	  w = sqrt (-2.0 * log (w) / w);

	  _smpl = x2 * w;
	  _havesmpl = true;

	  return x1 * w;
  }
}

complex<double> Random::complexRandn()
{
  double x1, x2, w, y1, y2;

  do
    {
      x1 = 2.0 * (rand_r (&_seed) / (double) RAND_MAX) - 1.0;
      x2 = 2.0 * (rand_r (&_seed) / (double) RAND_MAX) - 1.0;
      w = x1 * x1 + x2 * x2;
    }
  while (w >= 1.0);

  w = sqrt (-2.0 * log (w) / w);
  y1 = x1 * w;
  y2 = x2 * w;

  return complex<double> (y1, y2);
}


