#ifndef __SVTXCLUSTER_H__
#define __SVTXCLUSTER_H__

#include <phool/PHObject.h>
#include <vector>
#include <set>
#include <iostream>

class SvtxCluster : public PHObject {

public:
  
  typedef std::set<unsigned int>::const_iterator ConstHitIter;
  typedef std::set<unsigned int>::iterator       HitIter; 
  
  SvtxCluster();
  SvtxCluster(const SvtxCluster &clus);
  SvtxCluster& operator=(const SvtxCluster &clus);
  virtual ~SvtxCluster();

  // PHObject virtual overloads
  
  void         identify(std::ostream& os = std::cout) const;
  void         Reset();
  int          IsValid() const;

  // cluster info
  
  unsigned int get_id() const                        {return _id;}
  void         set_id(unsigned int id)               {_id = id;}
  
  unsigned int get_layer() const                     {return _layer;}
  void         set_layer(unsigned int layer)         {_layer = layer;}

  float        get_x() const                         {return _pos[0];}
  void         set_x(float x)                        {_pos[0] = x;}
  
  float        get_y() const                         {return _pos[1];}
  void         set_y(float y)                        {_pos[1] = y;}

  float        get_z() const                         {return _pos[2];}
  void         set_z(float z)                        {_pos[2] = z;}
  
  float        get_position(int coor) const          {return _pos[coor];}
  void         set_position(int coor, float xi)      {_pos[coor] = xi;}

  float        get_e() const                         {return _e;}
  void         set_e(float e)                        {_e = e;}

  unsigned int get_adc() const                       {return _adc;}
  void         set_adc(unsigned int adc)             {_adc = adc;}
  
  float        get_size(int i, int j) const;         //< get cluster dimension covar
  void         set_size(int i, int j, float value);  //< set cluster dimension covar

  float        get_error(int i, int j) const;        //< get cluster error covar
  void         set_error(int i, int j, float value); //< set cluster error covar

  //
  // clustered hit ids methods
  //
  void         clear_hits()                          {_hit_ids.clear();}
  bool         empty_hits()                          {return _hit_ids.empty();}
  size_t       size_hits()                           {return _hit_ids.size();}
  void         insert_hit(unsigned int hit_id)       {_hit_ids.insert(hit_id);}
  size_t       erase_hit(unsigned int hit_id)        {return _hit_ids.erase(hit_id);}
  ConstHitIter begin_hits() const                    {return _hit_ids.begin();}
  ConstHitIter find_hit(unsigned int hitid) const    {return _hit_ids.find(hitid);}
  ConstHitIter end_hits() const                      {return _hit_ids.end();}
  HitIter      begin_hits()                          {return _hit_ids.begin();}
  HitIter      find_hit(unsigned int hitid)          {return _hit_ids.find(hitid);}
  HitIter      end_hits()                            {return _hit_ids.end();}
  
  // deprecated interface
  
  float        get_phi_size() const;
  float        get_z_size() const;
  
private:
  
  unsigned int _id;                       //< unique identifier within container
  unsigned int _layer;                    //< detector layer id
  float _pos[3];                          //< mean position x,y,z
  float _e;                               //< cluster energy
  unsigned int _adc;                      //< cluster sum adc
  std::vector<std::vector<float> > _size; //< size covariance matrix (+/- cm^2)
  std::vector<std::vector<float> > _err;  //< error covariance matrix (+/- cm^2)
  std::set<unsigned int> _hit_ids;        //< list of cell hit ids 
  
  ClassDef(SvtxCluster, 1);
};

#endif
