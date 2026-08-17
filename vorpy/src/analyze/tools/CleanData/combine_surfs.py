import tkinter as tk
from tkinter import filedialog
import os


def read_off(file_path):
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.lstrip().startswith('#')]

    if not lines or lines[0] != 'OFF':
        raise ValueError(f"{file_path} is not a valid OFF file.")

    counts = lines[1].split()
    num_vertices = int(counts[0])
    num_faces = int(counts[1])

    vertices = lines[2:2 + num_vertices]
    faces = lines[2 + num_vertices:2 + num_vertices + num_faces]

    return vertices, faces


def combine_off_files(input_files, output_file):
    all_vertices = []
    all_faces = []

    vertex_offset = 0

    for file_path in input_files:
        vertices, faces = read_off(file_path)

        all_vertices.extend(vertices)

        for face in faces:
            parts = face.split()
            num_face_vertices = int(parts[0])

            indices = [int(i) + vertex_offset for i in parts[1:1 + num_face_vertices]]
            extra_data = parts[1 + num_face_vertices:]

            new_face = [str(num_face_vertices)] + [str(i) for i in indices] + extra_data
            all_faces.append(" ".join(new_face))

        vertex_offset += len(vertices)

        print(f"Added: {file_path}")
        print(f"    Vertices: {len(vertices)}")
        print(f"    Faces:    {len(faces)}")

    with open(output_file, 'w') as f:
        f.write("OFF\n")
        f.write(f"{len(all_vertices)} {len(all_faces)} 0\n")

        for vertex in all_vertices:
            f.write(vertex + "\n")

        for face in all_faces:
            f.write(face + "\n")

    print("\nCombined OFF file:")
    print(f"    Files:    {len(input_files)}")
    print(f"    Vertices: {len(all_vertices)}")
    print(f"    Faces:    {len(all_faces)}")
    print(f"    Output:   {output_file}")


def main():
    root = tk.Tk()
    root.withdraw()
    input_files = []
    while True:
        input_file = filedialog.askopenfilename(
            title="Select OFF file to combine",
            filetypes=[("OFF file", "*.off")]
        )
        if os.path.exists(input_file):
            input_files.append(input_file)
            print(input_file)
        else:
            break


    if not input_files:
        return

    output_file = filedialog.asksaveasfilename(
        title="Save combined OFF file",
        defaultextension=".off",
        filetypes=[("OFF files", "*.off"), ("All files", "*.*")]
    )

    if not output_file:
        return

    combine_off_files(input_files, output_file)


if __name__ == "__main__":
    main()