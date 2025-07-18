# FSSP_CA
This project is a graphical simulator for the classic computer science problem, the **Firing Squad Synchronization Problem (FSSP)**, implemented in Python. It provides a robust and user-friendly interface to visualize and analyze the behavior of cellular automata designed to solve this synchronization challenge.


The application is built using **Tkinter** for the graphical user interface (GUI), **Matplotlib** for rendering high-quality spacetime diagrams, and **NumPy** for efficient data manipulation.

---

## 🎯 About the Project

The Firing Squad Synchronization Problem asks for the design of a one-dimensional cellular automaton that can synchronize all its cells. Starting from a quiescent state, and triggered by a single "general" at one end, all cells must enter a special "fire" state at the exact same time for the first time.

This implementation successfully simulates this process using a time-optimal algorithm based on the principles developed by **Yves Mazoyer**.

### Key Features

*   **Optimal Time Solution:** Implements a sophisticated 8-state algorithm that achieves the optimal firing time of **T = 2n - 2** for a line of *n* cells.
*   **Interactive GUI:** A clean and intuitive interface built with Tkinter allows for easy control over the simulation.
*   **Dynamic Visualization:** Renders the evolution of the cellular automaton in real-time as a spacetime diagram using Matplotlib.
*   **Full Simulation Control:** Users can **Run**, **Stop**, and **Reset** simulations.
*   **Parameter Configuration:**
    *   Set the number of cells (**N**).
    *   Adjust the simulation speed via a slider.
    *   Configure the number of simulation steps per visual update for faster execution with large N.
*   **Advanced Viewing Tools:**
    *   **Zoom In / Zoom Out** functionality to inspect details of the spacetime diagram.
    *   **Smooth Scrolling** with the mouse wheel directly over the diagram.
*   **Data & Image Export:**
    *   **Save as PNG:** Export the final spacetime diagram as a high-resolution PNG image, perfect for reports and presentations.
    *   **Save as CSV:** Export the entire simulation history (the state of each cell at each time step) to a CSV file for further analysis in tools like MATLAB or Excel.

---

## 🚀 Getting Started

### Prerequisites

The application is written in **Python 3.6+**. The following libraries are required:

*   **NumPy:** For efficient array operations.
*   **Matplotlib:** For plotting the spacetime diagrams.

### Installation

1.  Clone the repository to your local machine:
    ```bash
    git clone https://github.com/your-username/your-repository-name.git
    cd your-repository-name
    ```

2.  It is highly recommended to use a virtual environment:
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
    ```

3.  Install the required packages using pip:
    ```bash
    pip install numpy matplotlib
    ```

OR: You can use Anaconda or something similar as well

### Running the Application

To start the simulator, simply run the main Python script from the project's root directory:

```bash
python fssp_gui.py 
