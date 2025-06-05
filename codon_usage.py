import tkinter as tk
from tkinter import messagebox, scrolledtext
from Bio.Seq import Seq
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def calculate_gc_content(dna):
    g = dna.count('G')
    c = dna.count('C')
    try:
        return round(((g + c) / len(dna)) * 100, 2)
    except ZeroDivisionError:
        return 0

def analyze_sequence():
    dna = text_input.get('1.0', tk.END).strip().upper().replace('\n', '').replace(' ', '')
    if not dna:
        messagebox.showerror("Error", "Please enter a DNA sequence.")
        return
    if any(base not in 'ATGC' for base in dna):
        messagebox.showerror("Error", "DNA sequence can only contain A, T, G, C.")
        return
    if len(dna) % 3 != 0:
        messagebox.showwarning("Warning", "DNA length is not multiple of 3. Results may be partial.")

    gc = calculate_gc_content(dna)
    mrna = Seq(dna).transcribe()
    protein = mrna.translate(to_stop=True)

    codon_usage = defaultdict(int)
    for i in range(0, len(mrna) - 2, 3):
        codon = str(mrna[i:i+3])
        codon_usage[codon] += 1

    aa_freq = Counter(str(protein))

    # Display textual results
    result_text.delete('1.0', tk.END)
    result_text.insert(tk.END, f"GC Content: {gc}%\n\nCodon Usage:\n")
    for codon, count in sorted(codon_usage.items()):
        result_text.insert(tk.END, f"{codon}: {count}\n")
    result_text.insert(tk.END, "\nAmino Acid Frequency:\n")
    for aa, freq in sorted(aa_freq.items()):
        result_text.insert(tk.END, f"{aa}: {freq}\n")

    # Plot amino acid frequency graph
    plot_aa_freq(aa_freq)

def plot_aa_freq(freq_dict):
    # Clear previous plot if any
    for widget in plot_frame.winfo_children():
        widget.destroy()

    amino_acids = list(freq_dict.keys())
    counts = [freq_dict[aa] for aa in amino_acids]

    fig, ax = plt.subplots(figsize=(6,3), dpi=100)
    ax.bar(amino_acids, counts, color='skyblue')
    ax.set_title('Amino Acid Frequency')
    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('Count')
    ax.grid(axis='y')

    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.draw()
    canvas.get_tk_widget().pack()

# Setup tkinter window
root = tk.Tk()
root.title("Codon Usage Analyzer by Aditya atul")
root.geometry("600x600")

label_input = tk.Label(root, text="Enter DNA Sequence (A,T,G,C):", font=("Arial", 12))
label_input.pack(pady=5)

text_input = scrolledtext.ScrolledText(root, height=6, width=70, font=("Courier", 12))
text_input.pack(pady=5)

analyze_btn = tk.Button(root, text="Analyze Sequence", command=analyze_sequence,
                        font=("Arial", 12), bg="black",  activebackground="gray")
analyze_btn.pack(pady=10)

label_result = tk.Label(root, text="Results:", font=("Arial", 12))
label_result.pack(pady=5)

result_text = scrolledtext.ScrolledText(root, height=10, width=70, font=("Courier", 11), bg="black", fg="lime")
result_text.pack(pady=5)

plot_frame = tk.Frame(root)
plot_frame.pack(pady=10, fill='both', expand=True)

root.mainloop()