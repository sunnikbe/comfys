#include "Particle.hpp"

// Constructor
Particle::Particle(double q_in, double m_in, arma::vec r_in, arma::vec v_in)
{
  q_ = q_in;
  m_ = m_in;
  r_ = r_in;
  v_ = v_in;
}

// Methods returning member variable values
double Particle::charge()
{
  return q_;
}

double Particle::mass()
{
  return m_;
}

arma::vec Particle::position()
{
  return r_;
}

arma::vec Particle::velocity()
{
  return v_;
}

// Method returning info
std::string Particle::info()
{
  // Make r and v vec strings
  std::ostringstream s, ss;
  r_.st().raw_print(s);
  v_.st().raw_print(ss);

  std::string info_ = "Charge: " + std::to_string(q_) + "\n"
                    + "Mass: " + std::to_string(m_) + "\n"
                    + "Position: " + s.str()
                    + "Velocity: " + ss.str() + "\n";

  return info_;
}
