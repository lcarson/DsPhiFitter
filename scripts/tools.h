// Include files
#include <dirent.h>
#include <errno.h>
#include <iostream>
int getdir(std::string dir, std::vector<std::string> &files){
  DIR *dp;
  struct dirent *dirp;
  if((dp  = opendir(dir.c_str())) == NULL) {
    std::cout << "Error(" << errno << ") opening " << dir << std::endl;
    return errno;
    }
  while ((dirp = readdir(dp)) != NULL) {
    std::string file(dirp->d_name);
        files.push_back(file);
  }
  closedir(dp);
  std::cout <<"\n" << std::endl;
  return 0;
}

void split(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiters = " ")
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    
    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}