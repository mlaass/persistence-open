#ifndef JSON_CONF_HPP
#define JSON_CONF_HPP

#include "picojson.h"

#include <string>
#include <fstream>
#include <streambuf>

/*specializations*/
namespace detail {
    template <typename t>
    t get(std::string key, picojson::object root)
     {
	 return root[key].get<t>();
     }
    template <>
    int get(std::string key, picojson::object root)
     {
	 return static_cast<int>(root[key].get<double>());
     }

    template<typename t>
    void inject(std::string key, const t &value, picojson::object &root)
    {
	  root[key] = picojson::value(value);
    }
};


class JSONConfig
{
public:
    picojson::object root;
    void load(std::string fname){
	std::ifstream t(fname);
	std::string json((std::istreambuf_iterator<char>(t)),
                 std::istreambuf_iterator<char>());
	picojson::value v;
	std::string err = picojson::parse(v, json);
	if (! err.empty()) 
	    throw(std::runtime_error("JSON parse: " + err ));
	
	if (! v.is<picojson::object>()) 
	    throw(std::runtime_error( "JSON is not an object"));

	root = v.get<picojson::object>();
		  
		  
  
    };

    template<typename t>
    t get(std::string key)
     {
	 return detail::get<t>(key,root);
     }

    // string injection, general inject needs detail:: specializations above. todo. not neede now.
    template<typename t>
    void set(std::string key, t value)
      {
	  detail::inject(key,value,root);
      }
    std::string dumps()
     {
	 return picojson::value(root).serialize();
     }
    void save(std::string filename)
     {
	 std::ofstream ofs(filename);
	 ofs << dumps() << std::endl;
     }	    
    
    
};
#endif
