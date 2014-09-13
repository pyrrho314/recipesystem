
class NIRI(DataClassification):
    name="NIRI"
    usage = "Applies to any data from the NIRI instrument."
    parent = "GEMINI"

    requirement = PHU(INSTRUME='NIRI')

newtypes.append(NIRI())
