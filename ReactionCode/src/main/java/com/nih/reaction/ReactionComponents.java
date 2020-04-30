package com.nih.reaction;

import java.util.List;
import java.util.Set;

public class ReactionComponents {

	Set<String> conditions;
	Set<String> agents;
	Set<String> catalysts;
	Set<String> solvents;
	int step;
	double temperature;
	double time;
	double pressure;
	double ph;
	List<String> yield;
	
	public Set<String> getConditions() {
		return conditions;
	}
	public void setConditions(Set<String> conditions) {
		this.conditions = conditions;
	}
	public Set<String> getAgents() {
		return agents;
	}
	public void setAgents(Set<String> agents) {
		this.agents = agents;
	}
	public int getStep() {
		return step;
	}
	public void setStep(int step) {
		this.step = step;
	}
	public double getTemperature() {
		return temperature;
	}
	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}
	public double getTime() {
		return time;
	}
	public void setTime(double time) {
		this.time = time;
	}
	public double getPressure() {
		return pressure;
	}
	public void setPressure(double pressure) {
		this.pressure = pressure;
	}
	public double getPh() {
		return ph;
	}
	public void setPh(double ph) {
		this.ph = ph;
	}
	public List<String> getYield() {
		return yield;
	}
	public void setYield(List<String> d) {
		this.yield = d;
	}
	public Set<String> getCatalysts() {
		return catalysts;
	}
	public void setCatalysts(Set<String> catalysts) {
		this.catalysts = catalysts;
	}
	public Set<String> getSolvents() {
		return solvents;
	}
	public void setSolvents(Set<String> solvents) {
		this.solvents = solvents;
	}

}
