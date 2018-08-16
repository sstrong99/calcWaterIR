#ifndef MAP_H
#define MAP_H

//abstract base class
class Map {
public:
Map(const float &cut,const float &Ddip): cut(cut),Ddip(Ddip),cut2(cut*cut) {};
  virtual ~Map() {};

  virtual float w(const float &E) const = 0;
  virtual float x(const float &w) const = 0;
  virtual float mud(const float &E) const = 0;
  virtual float p(const float &w) const = 0;
  virtual float kintra(const float &E1, const float &E2, const float &x1, const float &x2, const float &p1, const float &p2) const = 0;

  inline float getcut2() const { return cut2; };  //in nm, convert to A0
  inline float getcut()  const { return cut; };
  inline float getDdip() const { return Ddip; };

protected:
  const float cut, Ddip, cut2;
};

//inputs are in nm and ps
class MapAS2008 : public Map {
public:
  //constructor
  MapAS2008() : Map(0.7831,0.058) {};
  ~MapAS2008() {};

  inline float w(const float &E) const
  { return 3762.0 - 5060.0*E - 86225.0*E*E; };
  inline float x(const float &w) const
  { return 0.1934 - w*1.75e-5; };
  inline float mud(const float &E) const
  { return 0.1333 + 14.17*E; }; //see Yang and Skinner PCCP 12 p982, 2010
  inline float p(const float &w) const
  { return 1.611 + w*5.893e-4; };
  inline float kintra(const float &E1, const float &E2, const float &x1, const float &x2, const float &p1, const float &p2) const
  { return (-1789.0 + 23852.0*(E1+E2))*x1*x2 - 1.966*p1*p2; };
  //TODO: 1.966 is the value from the paper, it should be 2.508 instead
};

class MapGruenbaum2013 : public Map {
public:
  //constructor
MapGruenbaum2013() : Map(0.7831,0.067) {};
  ~MapGruenbaum2013() {};

  inline float w(const float &E) const
  { return 3760.2 - 3541.7*E - 152677.0*E*E; };
  inline float x(const float &w) const
  { return 0.19285 - w*1.7261e-5; };
  inline float mud(const float &E) const
  { return 0.1646 + 11.39*E + 63.41*E*E; };
  inline float p(const float &w) const
  { return 1.6466 + w*5.7692e-4; };
  inline float kintra(const float &E1, const float &E2, const float &x1, const float &x2, const float &p1, const float &p2) const
  { return (-1361.0 + 27165.0*(E1+E2))*x1*x2 - 1.887*p1*p2; };
};
#endif
