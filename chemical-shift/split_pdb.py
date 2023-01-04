import sys
import os

pname = 'Abeta40'

fread = open("./processed.pdb")
i = 0

fwrite = open("./%s_tmp.pdb" % pname, "a")
for lines in fread.readlines():
    if lines.split()[0] == "ENDMDL":
        fwrite.close()
        output = './pdbs/%s/%s_%d.pdb' % (pname, pname, i)
        os.system('mv ./%s_tmp.pdb %s' % (pname, output))
        i += 1
        fwrite = open("./%s_tmp.pdb" % pname, "a")
    else:
        fwrite.write(lines)
fread.close()

frames = 500
for j in range(frames):
    prefix = '%s_%d' % (pname, j * 49 + 1)
    input_data = './pdbs/%s/%s_%d' % (pname, pname, j * 49 + 1)
    os.system('cp %s.pdb ./%s_tmp.pdb' % (input_data, pname))
    os.system('tleap -s -f tleap.in')
    os.system('pmemd -O -i mini.in -p tmp.top -c tmp.crd -o tmp.mini.out -r tmp.mini.rst')
    os.system('cpptraj tmp.top analysis.in')
    os.system('grep -v TER ./tmp.md.pdb > %s_p.pdb' % (input_data))

# os.system('/home/faculty/hfchen/software/SPARTA+/sparta+ -in %s_p.pdb -out tmp_ab40_vae/%s.tab -outS
# tmp_ab40_vae/%s_struct.tab \ 1>tmp_ab40_vae/tmp_SPARTA+ 2>&1' % (input_data, prefix, prefix))
