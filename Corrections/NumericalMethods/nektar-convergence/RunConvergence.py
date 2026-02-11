import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.cElementTree as ET

def Generate_mesh(Nel):

    x0 = -1
    x1 = 1
    elements = Nel
    dx = (x1 - x0) / elements

    nodes = elements + 1

    nektar = ET.Element("NEKTAR")
    geometry = ET.SubElement(nektar, "GEOMETRY", DIM="1" ,SPACE="1")
    vertex = ET.SubElement(geometry, "VERTEX")

    for i in range(nodes):
        vertex_id = ET.SubElement(vertex, "V", ID=f"{i}").text=f"{x0 + i*dx:4.4f} {0:4.4f} {0:4.4f}"

    element = ET.SubElement(geometry, "ELEMENT")
    for i in range(elements):
        element_id = ET.SubElement(element, "S", ID=f"{i}").text=f"{i:8d} {(i+1):8d}"

    composite = ET.SubElement(geometry, "COMPOSITE")

    # Fluid domain
    composite_id = ET.SubElement(composite, "C", ID=f"0").text=f"S[0-{elements-1}]"
    # Left BC
    composite_id = ET.SubElement(composite, "C", ID=f"1").text=f"V[0]"
    # Right BC
    composite_id = ET.SubElement(composite, "C", ID=f"2").text=f"V[{nodes -1}]"

    # Domain ID
    domain = ET.SubElement(geometry, "DOMAIN")
    domain_id = ET.SubElement(domain, "D", ID=f"0").text=f"C[0]"

    tree = ET.ElementTree(nektar)
    ET.indent(tree, space="\t", level=0)
    tree.write(f"HelmholtzMesh_Nel{Nel}.xml")

def generate_prefinement(pstart, pend, pstep):

    P_arr = np.arange(pstart, pend, pstep)
    input_session = "HelmholtzSession.xml"

    pref_ndof = []
    pref_l2err = []
    pref_linferr = []
    pref_h1err = []

    # Running p-refinement
    for P in P_arr:

        # Prepare Pth-order session files
        print(f"Running P: {P}")
        output_session = f"HelmholtzSession-P{P}.xml"
        sed_cmd = f"sed 's/NUMMODES=\"[0-9]*\"/NUMMODES=\"{P+1}\"/' {input_session} > {output_session}"
        os.system(sed_cmd)

        # Get degrees of freedom
        dof_cmd = f"FieldConvert -m dof HelmholtzMesh.xml {output_session} stdout"
        res = subprocess.check_output(dof_cmd, shell=True, text=True)
        pref_ndof.append(res.split('Total number of DOF:', maxsplit=1)[-1].split(maxsplit=1)[0])

        # Run ADRSolver
        adr_cmd = f"ADRSolver HelmholtzMesh.xml {output_session}"
        res = subprocess.check_output(adr_cmd, shell=True, text=True)
        # print(res)
        pref_l2err.append(float(res.split('L 2 error (variable u) : ', maxsplit=1)[-1].split(maxsplit=1)[0]))
        pref_linferr.append(float(res.split('L inf error (variable u) : ', maxsplit=1)[-1].split(maxsplit=1)[0]))
        pref_h1err.append(float(res.split('H 1 error (variable u) : ', maxsplit=1)[-1].split(maxsplit=1)[0]))

    np.save('pref.npy', np.array([pref_ndof, pref_l2err, pref_linferr, pref_h1err]))

# Running h-refinement
def generate_hrefinement(hstart, hend, hstep):

    nel_arr = np.arange(hstart, hend, hstep)
    href_ndof = []
    href_l2err = []
    href_linferr = []
    href_h1err = []

    for nel in nel_arr:

        # Generate mesh 
        Generate_mesh(nel)
        print(f"Running Nel: {nel}")
        # output_session = f"HelmholtzSession-P{P}.xml"
        # sed_cmd = f"sed 's/NUMMODES=\"[0-9]*\"/NUMMODES=\"{P+1}\"/' {input_session} > {output_session}"
        # os.system(sed_cmd)

        # Get degrees of freedom
        dof_cmd = f"FieldConvert -m dof HelmholtzMesh_Nel{nel}.xml HelmholtzSession-P1.xml stdout"
        res = subprocess.check_output(dof_cmd, shell=True, text=True)
        href_ndof.append(res.split('Total number of DOF:', maxsplit=1)[-1].split(maxsplit=1)[0])

        # Run ADRSolver
        adr_cmd = f"ADRSolver HelmholtzMesh_Nel{nel}.xml HelmholtzSession-P1.xml"
        res = subprocess.check_output(adr_cmd, shell=True, text=True)
        # print(res)
        href_l2err.append(float(res.split('L 2 error (variable u) : ', maxsplit=1)[-1].split(maxsplit=1)[0]))
        href_linferr.append(float(res.split('L inf error (variable u) : ', maxsplit=1)[-1].split(maxsplit=1)[0]))
        href_h1err.append(float(res.split('H 1 error (variable u) : ', maxsplit=1)[-1].split(maxsplit=1)[0]))

    np.save('href.npy', np.array([href_ndof, href_l2err, href_linferr, href_h1err]))

if __name__ == '__main__':

    print("Hello world")

    # generate_prefinement(1, 17, 2)
    # generate_hrefinement(2, 120, 4)
    
    pref = np.load('pref.npy').astype(float)
    href = np.load('href.npy').astype(float)
    fig, ax = plt.subplots(1,2)

    ax[0].loglog(pref[0], pref[1], '-o', label=r'$p$-type refinement')
    ax[0].loglog(href[0], href[1], '-o', label=r'$h$-type refinement')
    ax[0].set_ylim([10**(-12), 10])
    ax[0].grid()
    ax[0].legend()
    ax[0].set_title(r'$(a)$')
    ax[0].set_xlabel('$N_{dof}$')
    ax[0].set_ylabel('$||\epsilon||_E$')

    ax[1].semilogy(pref[0], pref[3], '-o', label=r'$p$-type refinement')
    ax[1].semilogy(href[0], href[3], '-o', label=r'$h$-type refinement')
    ax[1].set_ylim([10**(-12), 10])
    # ax[1].set_xlim([0, 100])
    ax[1].grid()
    ax[1].legend()
    ax[1].set_title(r'$(b)$')
    ax[1].set_xlabel('$N_{dof}$')
    ax[1].set_ylabel('$||\epsilon||_E$')
    
    fig.set_size_inches(10,3)
    fig.subplots_adjust(wspace=0.3)
    fig.savefig('hp-convergence.pdf', bbox_inches='tight')
