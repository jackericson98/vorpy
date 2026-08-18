import tkinter as tk
from tkinter import filedialog, colorchooser


def recolor_off(input_file, output_file, rgb):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Remove blank/comment lines only for parsing the header
    data_indices = [i for i, line in enumerate(lines) if line.strip() and not line.lstrip().startswith('#')]

    if not data_indices or lines[data_indices[0]].strip() != 'OFF':
        raise ValueError("Selected file does not appear to be a valid OFF file.")

    counts_idx = data_indices[1]
    counts = lines[counts_idx].split()
    num_vertices = int(counts[0])
    num_faces = int(counts[1])

    # Find the vertex and face lines while allowing blank lines/comments
    geometry_indices = data_indices[2:]
    vertex_indices = geometry_indices[:num_vertices]
    face_indices = geometry_indices[num_vertices:num_vertices + num_faces]

    r, g, b = rgb
    rgba = f"{r:.4f} {g:.4f} {b:.4f} 1.0"

    for idx in face_indices:
        parts = lines[idx].split()

        if not parts:
            continue

        n = int(parts[0])
        face_data = parts[:n + 1]

        lines[idx] = " ".join(face_data) + f" {rgba}\n"

    with open(output_file, 'w') as f:
        f.writelines(lines)


def main():
    root = tk.Tk()
    root.withdraw()

    input_file = filedialog.askopenfilename(
        title="Select OFF file",
        filetypes=[("OFF files", "*.off"), ("All files", "*.*")]
    )

    if not input_file:
        return

    color = colorchooser.askcolor(title="Choose face color")

    if color[0] is None:
        return

    rgb = tuple(x / 255.0 for x in color[0])

    output_file = filedialog.asksaveasfilename(
        title="Save recolored OFF file",
        defaultextension=".off",
        filetypes=[("OFF files", "*.off"), ("All files", "*.*")]
    )

    if not output_file:
        return

    recolor_off(input_file, output_file, rgb)

    print(f"Input:  {input_file}")
    print(f"Color:  RGB{rgb}")
    print(f"Output: {output_file}")


if __name__ == "__main__":
    main()
