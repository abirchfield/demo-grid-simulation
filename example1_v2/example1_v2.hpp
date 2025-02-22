#include<stdio.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<time.h>

#include <sundials/sundials_core.h>
#include <idas/idas.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunlinsol/sunlinsol_dense.h>

#include "SystemModel.hpp"
#include "Bus.hpp"
#include "Bus.cpp"
#include "BusInfinite.hpp"
#include "BusInfinite.cpp"
#include "Branch.hpp"
#include "Branch.cpp"
#include "BusFault.hpp"
#include "GENROU.hpp"

#include "Ida.hpp"
#include "Ida.cpp"

int example1_v2();
