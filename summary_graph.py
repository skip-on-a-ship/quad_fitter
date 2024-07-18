import matplotlib.pyplot as plt

simulated_xs = [0,350,0,0,350,150,-350,0,-350,-350,350,-100]
simulated_ys = [0,0,350,0,350,0,0,-350,-350,350,-350,250]

reconstructed_xs = [-9.7,318.4,3.8,9.7,327.2,124.5,-305.4,-3.3,-342.9,-348.8,340.2,-75]
reconstructed_ys = [-13.1,-5.4,330.1,-4.8,331.2,0.7,2.4,-316.6,-343.6,350,-336,205.6]

plt.scatter(simulated_xs, simulated_ys, c="blue", label = "Simulated Locations", s=5)
plt.scatter(reconstructed_xs, reconstructed_ys, c="red", label = "Reconstructed Locations", s=5)

for i in range(len(simulated_xs)):
    plt.plot([simulated_xs[i], reconstructed_xs[i]],[simulated_ys[i], reconstructed_ys[i]], linewidth=0.7, color="black")

plt.legend()
plt.title("Simulated and Reconstructed Positions Around the Detector")
plt.xlim(-500,500)
plt.ylim(-500,500)
plt.savefig("graphs/summary_plot.pdf")
