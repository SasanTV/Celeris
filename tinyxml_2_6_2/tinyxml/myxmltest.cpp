#include "tinyxml.h"

#include <iostream>
#include <string>
#include <sstream>

using namespace std;

bool readInputCML(const char* pFilename)
{

	TiXmlDocument doc(pFilename);
	bool loadOkay = doc.LoadFile();
	if (loadOkay)
	{
		printf("\n%s:\n", pFilename);
		TiXmlElement* root = doc.FirstChildElement();
		if(root == NULL)
		{
			cerr << "Failed to load file: No root element."<< endl;
			doc.Clear();
			return 1;
		}
		for(TiXmlElement* elem = root->FirstChildElement(); elem != NULL; elem = elem->NextSiblingElement())
		{
			string elemName = elem->Value();
			const char* attr;
			if(elemName == "fieldDimensions")
			{
				cout << elemName << endl;
				float length,width;
				elem->QueryFloatAttribute("length", &length);
				elem->QueryFloatAttribute("width", &width);
				
				cout << length << " by " << width << endl;
			}
			 else if(elemName == "gridSize")
			 {

				cout << elemName << endl;
				float nx,ny;
				elem->QueryFloatAttribute("nx", &nx);
				elem->QueryFloatAttribute("ny", &ny);
				cout << nx << " by " << ny << endl;

			}
			else if(elemName == "bathymetryFilePath")
			{
				cout << elemName << endl;
				TiXmlText *t = elem->FirstChild()->ToText();
				
				cout << t->Value() << endl;

			}
			else if(elemName == "hotStartFilePath")
			{
				cout << elemName << endl;
				TiXmlText *t = elem->FirstChild()->ToText();
				
				cout << t->Value() << endl;

			}
		}

	}
	else
	{
		printf(doc.ErrorDesc(), pFilename);
	}
}

int main(void)
{
	readInputCML("example1.xml");
	string a;
	cin >> a;
	return 0;
}