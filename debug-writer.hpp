#ifndef DEBUGWRITER_HPP
#define DEBUGWRITER_HPP

namespace debug{
#ifndef DEBUG
    template<typename type> void write(std::string realm, type what)
    { // do nothing for debugging now
    }
#else
    template<typename type> void write(std::string realm, type what)
    {
	throw(std::runtime_error(realm + ": Not implemented"));
    }
    template<> void write(std::string realm, persistence::Curve what)
    {
	std::ofstream ofs("debug-"+realm+".csv");
	ofs << "val  index  x y" << std::endl; 
	for (const auto &c :what)
	{
	    ofs << c.val << " " << c.index << " " << c.vertex[0] << " " << c.vertex[1] << std::endl;
	}
    }

    template<> void write (std::string realm, persistence::Result result)
    {
/*
  struct Result{
    std::vector<Bar> bars;
    std::vector<Component> comps;
    std::vector<int> used;
    Curve pruned;
  };

*/
	
	{// write pruned curve
	    write(realm+"-pruned",result.pruned); 
	} 
	    
    }

    template<> void write(std::string realm, persistence::Extrema what)
    {
	{// write maxima
	    std::ofstream ofs("debug-"+realm+"-max.csv");
	    ofs << "index" << std::endl;
	    for (const auto &i:what.max) ofs << i << std::endl;
	    ofs.close();
	}
	{//write minima
	    std::ofstream ofs("debug-"+realm+"-min.csv");
	    ofs << "index" << std::endl;
	    for (const auto &i:what.min) ofs << i << std::endl;
	    ofs.close();
	}
	
	
    }
#endif // DEBUG
}
#endif
