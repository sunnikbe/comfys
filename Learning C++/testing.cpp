// declaration of variables
int i;
double d;
std::string s;

// assigning values to variables
i = 10;
d = 3.14;
s = "Hei på du";

// doing both the declaration and the assignment

int i = 10;
double d = 3.14;
std::string s = "Hei på du";

// declaring multiple variables
double x, y = 1.0, z;
// only y is assigned a value

// the const before a variable decleration makes it impossible to accidentaly
// change a variable. If yu try to change it, an error messange will occur when
// compiling

const double g = 9.81 // [m/pow(s,2)] power function is pow(Base,Power)

// a simple code
std::string s; // if this part is inside the if block, or other curly bracket
// {} blocks, the code won´t work.

if (true)
{
  s = "Hello!";
}

std::cout << s;
