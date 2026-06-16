# Network visualization tutorial

## Quick Overview

<p align="justify">RingNet is an interactive network visualization platform for multi modal biological data. Each node can display up to five aligned data layers. The frontend loads JSON network objects and provides filtering, exploration, pattern discovery, and export functions.
<p>

## 1. Network control panel
<img width="1046" height="65" alt="image" src="https://github.com/user-attachments/assets/76446a59-0bc9-4d36-9dde-fbbe0c3a1551" />

<p align="justify">
This section controls network selection, layout, degree filtering, search, and export functions. Users can select networks, adjust layout, define degree range filters, and export results in PNG, SVG, or JSON formats.
<p>

### Network Selection

<p align = "justify">Network Selection: Users select a network from a dropdown menu. Only one network is rendered at a time. This improves performance and clarity. Select the specific network or community you wish to explore. Note that **only one network** is displayed at a time to maintain clarity.
<p>
  
### Layout Adjustment
Arrange nodes on the canvas using various algorithms. Changing the layout **only affects node positions** and does not modify the underlying data.

| Layout Type | Best Use Case |
| :--- | :--- |
| **Loose Balanced / Soft Spread** | General network exploration |
| **Spread Force** | Highlighting clusters |
| **Component Grid** | Separating disconnected components |
| **Concentric / Circle** | Hierarchical or symmetric structures |
| **Grid** | Uniform distribution |

### Degree Range Filter
Filter nodes based on their connectivity.
> **Example:** Setting `Degree Range: 0 – 500` will only show nodes whose visible degree falls within this specific range.

### Export & Snapshot Management
<img width="1056" height="75" alt="image" src="https://github.com/user-attachments/assets/cda86c6f-7ec6-4d19-ad87-b49e994ee48a" />

Keep your work safe and ready for publication with these tools:

- **Fit:** Instantly fit the entire graph into the visible canvas.
- **Save PNG:** Export the current view as a standard image.
- **Save SVG:** Export as an editable vector image (Recommended for **publications**).
- **Save JSON:** Export the raw network data in JSON format.
- **Save Snapshot (Highly Recommended):** Save the entire visualization setting customized by the user in a **zip file**  and the user can upload the **zip file** to restore the view. This includes:
  - Node positions, colors, and filters.
  - Pattern Discovery and Similarity/Difference settings.
  - Highlighted nodes/edges and search states.

### Node Search & Interaction

- **Node Search:** Search for specific genes or node names. The matched node and its immediate neighbors will be highlighted, while others are dimmed.
- **Clear:** Removes search highlights. If **Pattern Discovery** is active, clearing the search returns you to the pattern state rather than resetting the entire graph.

## 2. Data Layer Configuration

<img width="904" height="68" alt="image" src="https://github.com/user-attachments/assets/a0014191-d35f-44b6-9d47-5a116b087242" />

This section defines how multi modal data is visualized using concentric rings within nodes.

### Color Palette:
We provide four color palettes including the style from the New England Journal of Medicine (NEJM), Nature Cancer, Lancet, and Science for users to choose (see Section 6 for details).

### Different data types:

| Layer | Common Biological Meaning |
| :--- | :--- |
| **Data1** | Continuous data type |
| **Data2** | Continuous data type |
| **Data3** | Discrete data type |
| **Data4** | Discrete data type  |
| **Group** | Stage, Grade, or Class labels |

> **Customization:** Each layer supports a three-color gradient: 
> `[ Negative Color ]` — `[ Neutral Color ]` — `[ Positive Color ]`

### Data Filtering & Normalization

Refine your view by selecting specific data layers and scales:

1. **Select Layer:** Choose between Data1, Data2, Data3, or Data4.
2. **Choose Scale:** 
   - `Raw`: Original values.
   - `Min-Max Norm`: Scaled between 0 and 1.
   - `Z-score`: Standardized relative to the mean.
3. **Value Range Filtering:**  Users can filter nodes based on value ranges for selected data layers.

## 3. Pattern Discovery Module

<img width="1844" height="152" alt="image" src="https://github.com/user-attachments/assets/587d11bc-a901-46af-abff-72c2a25ac81f" />

<p align="justify">
Pattern Discovery identifies anti correlated multi omics patterns such as high expression with low methylation. Thresholds define high and low states, while intensity is computed using z score based measures
<p>
  
### How it works:
Pattern Discovery typically compares **Data1** and **Data2** using Z-score thresholds:
- **Data1 High/Low:** Define thresholds for high/low states.
- **Data2 High/Low:** Define thresholds for high/low states.
- **Highlighting:** Nodes matching opposite directions trend (High/Low or Low/High) are automatically highlighted.

### Analysis Modes (Group-aware vs. Group-free)

<img width="1842" height="63" alt="image" src="https://github.com/user-attachments/assets/eaea3244-9901-4ee3-bd93-2438c2e921f9" />

Depending on whether you have stage or group information, the viewer adjusts its metrics:

### **Group-aware Mode** (With Stage Data)
- **TrendScore:** Measures if the pattern follows a group-level trend.
- **GroupScore:** Measures the strength of group-aware evidence.
- *Best for: Stage, Class, or Condition analysis.*

### **Group-free Mode** (No Stage Data)
- **GlobalPatternRate:** Measures how many samples support the pattern globally.
- **GroupScore:** Measures global evidence without using labels.


## 4. Difference and Similarity Displaying Module

<img width="441" height="60" alt="image" src="https://github.com/user-attachments/assets/2e9508b4-3550-46e9-9c3c-9d68bf42dc5d" />

This module selects representative samples based on multi-omics deviation. It supports stage aware sampling and robust node filtering using Data1 and Data2 thresholds.

- **Difference Displaying:** Selects samples with the most multi-omics variation.
- **Similarity Displaying:** Selects samples that are most coherent.
- **Samples per Group:** Controls the number of selected samples (Default: `8`).


## 5. Group and Edge Visualization Controls

<img width="904" height="92" alt="image" src="https://github.com/user-attachments/assets/9a3969a1-44b6-4c2b-a2b2-705c1fbd1fed" />

Groups and Edges are color coded, and users can define their colors. Users can also filter edges by weight and interaction type.
Fine-tune the connections between nodes:

- **Edge Colors:** Map weights to a color gradient (Negative/Neutral/Positive).
- **Show Weights:** Toggle numeric labels on edges.
- **Edge Data Norm:** Choose between `Raw`, `Min-Max`, or `Z-score`.
- **Color by `interact_id`:** Color edges by interaction type.

## 6. Color Palette 

<img width="1038" height="101" alt="image" src="https://github.com/user-attachments/assets/f7471bb6-1fa0-4ceb-952a-8790943b682c" />

Only one palette can be selected at a time. The palette checkboxes are mutually exclusive. Selecting one palette automatically deselects the previously selected palette. Unchecking the currently selected palette restores the original RingNet colors. Manually changing any Data or Edge color automatically clears the selected palette. 

The RingNet palette uses eight-digit hexadecimal color values. For example, in the format #RRGGBBAA, RR represents the red channel, GG represents the green channel, BB represents the blue channel, and AA represents the alpha or transparency channel. For example, in the value #0072B5FF, 00 is the red value, 72 is the green value, B5 is the blue value, and FF indicates a solid color. Because the browser color picker uses six-digit hexadecimal values, RingNet displays the RGB portion in the color input while preserving the intended fully opaque appearance.

The palette controls are displayed above the Data 1–Data 4 color selectors.

### Available options:

- **Default color palette**

<img width="904" height="68" alt="image" src="https://github.com/user-attachments/assets/a0014191-d35f-44b6-9d47-5a116b087242" />

When no journal palette is selected, RingNet uses its original default colors including blue, white, and red. The original Data 1–Data 4 color selectors and Edge color selectors remain editable. To restore the original colors: Uncheck the active journal palette. RingNet automatically returns to the default color configuration.

- **The New England Journal of Medicine (NEJM)**

<img width="687" height="70" alt="image" src="https://github.com/user-attachments/assets/4a5a83c0-518c-4cbb-8189-4b010ab2e231" />

Color configuration:

Data 1: [#0072B5FF, #FFFFFFFF, #BC3C29FF]

Data 2: [#20854EFF, #FFFFFFFF, #E18727FF]

Data 3: [#FFDC91FF, #FFFFFFFF, #7876B1FF]

Data 4: [#FFFFFFFF, #FFFFFFFF, #6F99ADFF]

Edge: [#7FBC41FF, #F7F7F7FF, #DE77AEFF]

- **Nature Cancer**

<img width="687" height="64" alt="image" src="https://github.com/user-attachments/assets/3bdf568e-8b1c-4e22-be34-61956a1500d9" />

Color configuration:

Data 1: [#4DBBD5FF, #FFFFFFFF, #E64B35FF]

Data 2: [#00A087FF, #FFFFFFFF, #7E6148FF]

Data 3: [#F39B7FFF, #FFFFFFFF, #8491B4FF]

Data 4: [#FFFFFFFF, #FFFFFFFF, #1A1A1AFF]

Edge: [#D6604DFF, #F7F7F7FF, #878787FF]

- **Lancet**

<img width="684" height="67" alt="image" src="https://github.com/user-attachments/assets/b754bc87-bdb1-4195-991e-ac886437b829" />

Color configuration:

Data 1: [#00468BFF, #FFFFFFFF, #ED0000FF]

Data 2: [#42B540FF, #FFFFFFFF, #FDAF91FF]

Data 3: [#0099B4FF, #FFFFFFFF, #925E9FFF]

Data 4: [#FFFFFFFF, #FFFFFFFF, #1B1919FF]

Edge: [#42B540FF, #F7F7F7FF, #AD002AFF]

- **Science**

<img width="687" height="65" alt="image" src="https://github.com/user-attachments/assets/a81c6d64-7be5-4d36-acc9-59a870407ecf" />

Data 1: [#3B4992FF, #FFFFFFFF, #EE0000FF]

Data 2: [#008B45FF, #FFFFFFFF, #A20056FF]

Data 3: [#B35806FF, #FFFFFFFF, #631879FF]

Data 4: [#FFFFFFFF, #FFFFFFFF, #1B1919FF]

Edge: [#008280FF, #F7F7F7FF, #BB0021FF]

### Directed and undirected networks

The palette behavior is identical in directed and undirected network viewers. The palette changes only the visual color mapping. It does not change:

- Node values
- Edge weights
- Network topology
- Direction information
- Community structure
- Data normalization
- Filtering thresholds
- Pattern Discovery results
- Difference displaying results
- Similarity displaying results

For directed networks, arrow direction and directed edge structure are preserved. For undirected networks, merged edge structure and undirected rendering are preserved.

### Saving palette settings

The currently selected palette is included in the viewer configuration. When using: `Save JSON` or `Save Snapshot`, RingNet stores the active palette selection and the current visualization settings together. When the configuration is restored, the corresponding palette is restored as well.

## 7. Viewer Types

| Directed Viewer | Undirected Viewer |
| :--- | :--- |
| Preserves edge direction (`A → B` ≠ `B → A`). | Ignores direction (`A → B` & `B → A` are merged). |
| Best for: Regulatory or Signaling networks. | Best for: Similarity or Co-expression networks. |
