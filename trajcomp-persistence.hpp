#ifndef TRAJCOMP_PERSISTENCE_HPP
#define TRAJCOMP_PERSISTENCE_HPP

/*
trajcomp cleaned implementation of persistence simplification and friends
*/

#include<trajcomp/trajcomp.hpp>

namespace persistence{

    // TODO: remove all this helper functions and replace with typical generics like
    // glm
      typedef std::vector<double> vec2;
  typedef std::vector<vec2> Trajectory;
  
  double lengthSqr(const vec2 & v){
    return v[0] * v[0] + v[1] * v[1];
  }
  
  double length(const vec2 & v){
    return sqrt(lengthSqr(v));
  }
  
  double dot(const vec2 & v, const vec2 & w){
    return v[0] * w[0] + v[1] * w[1];
  }
  double det(const vec2 & v, const vec2 & w){
    return v[0] * w[1] + v[1] * w[0];
  }
  
  //positive if v is left of w negative else
  double left_right(const vec2 & v, const vec2 & w){
    return v[0] * w[1] - v[1] * w[0];
  }
  
  //1.0 if v is left of w -1 else
  double lr_sgn(const vec2 & v, const vec2 & w){
    return (left_right(v, w) < 0) ? -1.0: 1.0;
  }
  
  double to_deg(double rad){
    return (rad * 180.0) / M_PI;
  }

  vec2 mult(const vec2 & v, double x){
    return {v[0] * x, v[1] * x};
  }
  
  vec2 normalize(vec2  v){
    auto len = length(v);
    v[0] *= 1/len;
    v[1] *= 1/len;
    return v;
  }
  
  vec2 diff( const vec2 & v, const vec2 & w){
    return {v[0] - w[0], v[1] - w[1]};
  }
  
  double distance(const vec2 & v, const vec2 & w){
    return length(diff(v,w));
  }

//////// end sensible vector helpers
struct CurveElement{
    double val;
    int index;
    vec2 vertex;
};

typedef std::vector<CurveElement> Curve;

  Curve traj_to_curve(const std::vector<vec2> &t){
      Curve ret;
      for (auto i=0; i<t.size(); i++)
      {
	  if (i == 0 || i == t.size() -1)
	  {
	      ret.push_back({0,i,t[i]});
	      continue;
	  }
	  // inner point
	  
	  
	  auto v = t[i];
	  auto last = t[i-1];
	  auto v1 = normalize(diff(t[i-1],t[i]));
	  //auto next = t[i+1];
	  auto v2 = normalize(diff(t[i], t[i+1]));
	  auto x = dot(v1, v2);
	  auto y = det(v1,v2);
	  auto degree = to_deg(atan2(y,x));

	  if(std::isnan(degree))
	      continue; // should be safe, but i does not occur anymore as a curve-element
	  ret.push_back({degree, i, t[i]});	
	  
      }
    return ret;  
  };
   

/*

Finding Extrema

 */    


      struct Extrema{
    std::vector<int> min, max;
  };
  
  Extrema extrema(const Curve &curve){
    
    const int fl_un = 0, fl_min = 1, fl_max = 2;
    Extrema res;
    if (curve.size() <=2)
    {
      if(curve.size() == 0)
        throw(std::runtime_error("Trying to compute extrema from trajectory without points"));

      res = {{0}, {static_cast<int>(curve.size()-1)}}; // can be 0,0
      return res;
    }
    auto d_min = curve[0], d_max = curve[0];
    auto search = fl_un;
    auto i_un = 0;
    
    for(auto i = 1; i < curve.size(); ++i){
      auto x = curve[i];
      if(search == fl_un){
        if(x.val > d_min.val){
          search = fl_max;
          d_max = x;
          res.min.push_back(i_un);
        }
        if(x.val < d_min.val){
          search = fl_min;
          d_min = x;
          res.max.push_back(i_un);
        }
      } else if(search == fl_min){
        if(x.val > d_min.val){
          search = fl_max;
          d_max = x;
          res.min.push_back(i-1);
        }else{
          d_min = x;
        }
      } else if(search == fl_max){
        if(x.val < d_max.val){
          search = fl_min;
          d_min = x;
          res.max.push_back(i-1);
        }else{
          d_max = x;
        }
      }
    }
    // if res.min and res.max are empty
    if (res.min.size() == 0 || res.max.size() == 0)
    {
      res.min.push_back(0);
      res.max.push_back(curve.size()-1);
      //	throw(std::runtime_error("Unable to decide on end as we did not record any extrema yet"));
    }
    
    //add end, its either min or max
    if(res.min.back() > res.max.back()){
      res.max.push_back(curve.size()-1);
    }else{
      res.min.push_back(curve.size()-1);
    }
    return res;
  }

/*
   BARS and Pruning
 */

  struct Component{
    int left, right;
    int min, max;
    bool finished;
  };
  
  struct Bar{
    int min, max;
  };
  
  struct Result{
    std::vector<Bar> bars;
    std::vector<Component> comps;
    std::vector<int> used;
    Curve pruned;
  };
    
  Curve betaPruning(const std::vector<Bar> &bars, const Curve &curve, double beta = 0.01){
    std::set<CurveElement, std::function<bool (const CurveElement&, const CurveElement&)>> pruned([](const CurveElement& a, const CurveElement& b){
        return a.index < b.index;
      });
    
    for(auto b : bars){
      if(curve[b.max].val - curve[b.min].val > beta){
        pruned.insert(curve[b.min]);
        pruned.insert(curve[b.max]);
      }
      // else{
      //   std::cout <<"pruned: ["<<b.min<<", "<<b.max<<"]\n";
      // }
    }
    Curve pruned_curve;
    std::copy(pruned.begin(), pruned.end(), std::back_inserter(pruned_curve));
    return pruned_curve;
  }

/*
THE CORE ALGORITHM
 */

  Result persistence(const Curve &curve, double beta = 0.01, int debug_iter = -1){
    
    auto ex = extrema(curve);
    auto min = ex.min;
    auto max = ex.max;
    
    auto used = std::vector<int>(curve.size());
    std::fill(used.begin(), used.end(), -1);
    
    auto max_set = std::set<int>(max.begin(), max.end());
    
    auto comps = std::vector<Component>();
    for(auto m:min){
      comps.push_back({m,m,m,-1, false});
    }
    auto active_comps = std::vector<int>();
    for(int i = 0; i < comps.size(); ++i){
      active_comps.push_back(i);
    }
        
    //setup loop
    bool finished = false;
    int iter = -1;
    if(debug_iter < 0){
      debug_iter = curve.size()*(int)comps.size();
    }
    
    int max_iter = std::min((int)curve.size()*(int)comps.size(), debug_iter);
    
    while (! finished && ++iter < max_iter && active_comps.size() > 0){      
      //grow bars
      int i = iter%((int) active_comps.size());
      //std::cout<<"iter: "<<iter << " %% "<< i<<"/"<<active_comps.size() <<std::endl;
      auto current_i = active_comps[i];
      auto &c1 = comps[current_i];
      
      if(c1.left == 0 && c1.right == (int)curve.size()-1){        
        finished = true; 
        c1.finished = true;
        //REMOVE
        active_comps.erase(active_comps.begin()+i);
        
        if(max_set.count(c1.left) >0){
          c1.max = c1.left;
        } else if(max_set.count(c1.right) >0){
          c1.max = c1.right;
        }
        
///        std::cout<<"Finished!!!!!!!!!!!!!!!! at: "<< iter<<std::endl;
      }
      
      if(!c1.finished){
        auto l = c1.left - 1;
        auto r = c1.right + 1;
        
        //decide grow direction
        int x = 0;
        if((curve[l].val < curve[r].val || r == curve.size())  && l != -1  ){
          c1.left = l;
          x = l;
        }else{
          c1.right = r;
          x = r;
        }
        
        if(x== -1){
          x=0;
        }

        if(x == curve.size()){
          x = curve.size()-1;
        }
        
        //if is _max(x)
        if(max_set.count(x) > 0 ){

          c1.finished = true;
          c1.max = x;
          //REMOVE
          active_comps.erase(active_comps.begin()+i);

          //if max is used
          if(used[x] != -1){
            
            //merge
            auto & c2 = comps[used[x]];
            auto m1 = c1.min;
            auto m2 = c2.min;
            if (c2.min > c1.min){
              m1 = c2.min;
              m2 = c1.min;
            }
            
            c2.left = std::min(c1.left, c2.left);
            c2.right = std::max(c1.right, c2.right);
            c2.min = m2;
            c2.max = -1;
            c2.finished = false;
            //ADD
            active_comps.push_back(used[x]);
            
            c1.left = std::min(x, m1);
            c1.right = std::max(x, m1);
            c1.min = m1;
            //REDUNDANT see above
            c1.max = x;
            c1.finished = true;

            ////REMOVE
            //active_comps.erase(active_comps.begin()+i);
          }
          used[x] = current_i;
        }
      }
    }
    
    //make bars
    std::vector<Bar> bars;
    
    for(auto c : comps){
      if(c.max >=0){
        bars.push_back({c.min,c.max});
      }
    }

    auto pruned_result = betaPruning(bars, curve, beta);
    if (pruned_result.size() == 0)
    {
	    pruned_result = {curve[0], curve.back()};
    }else{
      if (pruned_result[0].index != 0)
      {
      // std::cout << "Patching start into result " << std::endl;
          pruned_result.insert(pruned_result.begin(), curve[0]); 
      }
      
      if (pruned_result.back().index != 0)
      {
      // std::cout << "Patching end into result " << std::endl;
          pruned_result.push_back(curve.back()); 
      }
    }
    
    // hack: enforce start and end to be part of the result.
    return {bars, comps, used, pruned_result};
  };


/*
PRUNING: Hierarchical first
*/
  Curve prune_curve_dist(const Curve &curve, double scale){
    Curve res;
    if (curve.size() < 2) // empty curve prunes to empty curve
       return res;
    res.reserve(curve.size());
    res.push_back(curve[0]);
    
    for(auto i = 1; i < curve.size()-1;  ++i){
      auto it = curve[i], it2 = curve[i+1];
      //check distance
      if(distance(it.vertex, it2.vertex) < scale){
        //pruning
        if(it.val < it2.val){
          //avoid double insertion
          if(it.index != res.back().index)
            res.push_back(it);
        }else {
          res.push_back(it2);
        }
      }else{
        //no pruning
        res.push_back(it);
      }
    }
    res.push_back(curve.back());
    return res;
  }


  Curve prune_sds(const Curve &curve, double eps){
    Curve res;
    if (curve.size() < 2) // empty curve prunes to empty curve
       return res;
    res.reserve(curve.size());
    res.push_back(curve[0]);
    trajcomp::default_segment_distance<vec2> sd;
    
    for(auto i = 0; i < curve.size()-2;  ++i){
	const auto &a = curve[i].vertex;
    	const auto &b = curve[i+1].vertex;
	const auto &c = curve[i+2].vertex;
	double d = sd(a,c,b);
	if ( d > eps)
	{
	    res.push_back(curve[i+1]);
	}
    }
    res.push_back(curve.back());
    return res;
  }

    
  Curve persistenceMultiRes(const Curve &curve,  double beta, int levels){
	  auto _curve = curve;
	  for(int i = 0; i < levels; ++i){
	      auto p_result = persistence(_curve, beta, -1);
	      _curve = prune_curve_dist(p_result.pruned,  pow(2,i));
	  }
	  return _curve;
  }


    
}; //namespace


#endif
