#include "writeSurfFields.H"
#include "OFstream.H"
#include "floatScalarSB.H"
#include "vtkWriteFieldOps.H"
#include "emptyFvsPatchFields.H"
#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::writeSurfFields
(
    const bool binary,
    const vtkMesh& vMesh,
    const fileName& fileName,
    const UPtrList<const surfaceVectorField>& surfVectorFields
)
{
    const fvMesh& mesh = vMesh.mesh();

    std::ofstream str(fileName.c_str());

    vtkWriteOps::writeHeader
    (
        str,
        binary,
        "surfaceFields"
    );

    str << "DATASET POLYDATA" << std::endl;

    const pointField& fc = mesh.faceCentres();

    str << "POINTS " << mesh.nFaces() << " float" << std::endl;

    DynamicList<floatScalar> pField(3*mesh.nFaces());

    for (label facei = 0; facei < mesh.nFaces(); facei++)
    {
        vtkWriteOps::insert(fc[facei], pField);
    }

    vtkWriteOps::write(str, binary, pField);

    str << "POINT_DATA " << mesh.nFaces() << std::endl
        << "FIELD attributes " << surfVectorFields.size() << std::endl;

    // surfVectorFields
    forAll(surfVectorFields, fieldi)
    {
        const surfaceVectorField& svf = surfVectorFields[fieldi];

        str << svf.name() << " 3 "
            << mesh.nFaces() << " float" << std::endl;

        DynamicList<floatScalar> fField(3*mesh.nFaces());

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            vtkWriteOps::insert(svf[facei], fField);
        }

        forAll(svf.boundaryField(), patchi)
        {
            const fvsPatchVectorField& pf = svf.boundaryField()[patchi];

            const fvPatch& pp = mesh.boundary()[patchi];

            if (isA<emptyFvsPatchVectorField>(pf))
            {
                // Note: loop over polypatch size, not fvpatch size.
                forAll(pp.patch(), i)
                {
                    vtkWriteOps::insert(vector::zero, fField);
                }
            }
            else
            {
                forAll(pf, i)
                {
                    vtkWriteOps::insert(pf[i], fField);
                }
            }
        }

        vtkWriteOps::write(str, binary, fField);
    }
}


// ************************************************************************* //
