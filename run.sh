echo "Hello :)"

echo "Download data files"
wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m11.txt
wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m12.txt
wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m13.txt
wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m14.txt
wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m15.txt


echo "a" 
python3 a.py

echo "b" 
python3 b_ref.py

echo "c" 
python3 c.py

echo "d" 
python3 d.py

echo "e" 
python3 e.py


echo "Generating the pdf"
pdflatex template.tex
bibtex template.aux
pdflatex template.tex
pdflatex template.tex
