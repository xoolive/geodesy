#ifndef LOGGER_H
#define LOGGER_H

#include <sstream>
#include <iostream>
#include <cmath>
#include <cfloat>

class Logger
{
  std::ostringstream out;
  std::string file;
  int nbError;

public:
  Logger(std::string file) : file(file), nbError(0) {}

  void message (int line, std::string message)
  {
    out << line << ":\t" << message << std::endl;
  }

  void testint (int line, int a, int b, std::string message)
  {
    out << line << ":\t" <<
      "a=" << a <<"; b=" << b  << "\t" << message << std::endl;
    if (a != b)
    {
      out << "ERROR ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ERROR\n";
      nbError += 1;
    }
  }

  void testfloat (int line, float a, float b, std::string message,
                  float precision = FLT_EPSILON)
  {
    out << line << ":\t" <<
      "a=" << a <<"; b=" << b  << "\t" << fabs(a - b) << "\t" <<
      message << std::endl;
    if (fabs(a - b) > precision)
    {
      out << "ERROR ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ERROR\n";
      nbError += 1;
    }
  }

  void testdouble (int line, double a, double b, std::string message,
                  double precision = DBL_EPSILON)
  {
    out << line << ":\t" <<
      "a=" << a <<"; b=" << b  << "\t" << fabs(a - b) << "\t" <<
      message << std::endl;
    if (fabs(a - b) > precision)
    {
      out << "ERROR ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ERROR\n";
      nbError += 1;
    }
  }

  void testhex (int line, int a, int b, std::string message)
  {
    out << line << ":\t" << std::hex <<
      "a=0x" << a <<"; b=0x" << b << std::dec <<
      "\t" << message << std::endl;
    if (a != b)
    {
      out << "ERROR ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ERROR\n";
      nbError += 1;
    }
  }

  int reportexit() {
    std::cerr << "############################### " << nbError << " error(s) in: \n" << file << std::endl;
    std::cerr << out.str();
    if (0 == nbError) return EXIT_SUCCESS;
    std::cerr << "############################### " << nbError << " error(s) in: \n" << file << std::endl;
    return EXIT_FAILURE;
  }

};

#endif // LOGGER_H

