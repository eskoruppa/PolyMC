#include "Argparse.h"


std::string parse_arg(const std::string& default_val,const std::string &keyword, const std::vector<std::string> & argv) {
    int argc = argv.size();
    std::string argstr;
    for (int i=1;i<argc;i++) {
		argstr = argv[i];
		if (argstr==keyword && argc > i) {
            argstr = argv[i+1];
            return argstr;
		}
    }
    return default_val;
}

bool parse_flag(const std::string &keyword, const std::vector<std::string> & argv) {
    int argc = argv.size();
    std::string argstr;
    for (int i=1;i<argc;i++) {
		argstr = argv[i];
		if (argstr==keyword) {
            return true;
		}
    }
    return false;
}

bool parse_flag(bool default_val, const std::string &keyword, const std::vector<std::string> & argv) {
    int argc = argv.size();
    std::string argstr;
    for (int i=1;i<argc;i++) {
		argstr = argv[i];
		if (argstr==keyword) {
            return !default_val;
		}
    }
    return default_val;
}

bool parse_flag_atpos(const std::string &keyword, int flag_pos, const std::vector<std::string> & argv) {
    int argc = argv.size();
    std::string argstr;
    if (argc >= flag_pos+1) {
        argstr = argv[flag_pos];
		if (argstr==keyword) {
            return true;
		}
    }
    return false;
}

std::string parse_str_atpos(int pos, const std::vector<std::string> & argv) {
    int argc = argv.size();
    std::string argstr;
    if (argc >= pos+1) {
        argstr = argv[pos];
        return argstr;
    }
    return "";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//std::string parse_arg(const std::string& default_val,const std::string &keyword, int argc, const char ** argv) {
//    std::string argstr;
//    for (int i=1;i<argc;i++) {
//		argstr = argv[i];
//		if (argstr==keyword && argc > i) {
//            argstr = argv[i+1];
//            return argstr;
//		}
//    }
//    return default_val;
//}
//
//bool parse_flag(const std::string &keyword, int argc, const char ** argv) {
//    std::string argstr;
//    for (int i=1;i<argc;i++) {
//		argstr = argv[i];
//		if (argstr==keyword) {
//            return true;
//		}
//    }
//    return false;
//}
//
//bool parse_flag(bool default_val, const std::string &keyword, int argc, const char ** argv) {
//    std::string argstr;
//    for (int i=1;i<argc;i++) {
//		argstr = argv[i];
//		if (argstr==keyword) {
//            return !default_val;
//		}
//    }
//    return default_val;
//}
//
//bool parse_flag_atpos(const std::string &keyword, int flag_pos , int argc, const char ** argv) {
//    std::string argstr;
//    if (argc >= flag_pos+1) {
//        argstr = argv[flag_pos];
//		if (argstr==keyword) {
//            return true;
//		}
//    }
//    return false;
//}
//
//std::string parse_str_atpos(int pos , int argc, const char ** argv) {
//    std::string argstr;
//    if (argc >= pos+1) {
//        argstr = argv[pos];
//        return argstr;
//    }
//    return "";
//}
